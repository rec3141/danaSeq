<script>
  let { status = null } = $props();

  // Timeline toggle
  let showTimeline = $state(false);

  // Time since last update
  let now = $state(Date.now());
  $effect(() => {
    const timer = setInterval(() => { now = Date.now(); }, 5000);
    return () => clearInterval(timer);
  });

  let updatedAgo = $derived.by(() => {
    if (!status?.timestamp) return '';
    const diff = Math.floor((now - new Date(status.timestamp).getTime()) / 1000);
    if (diff < 10) return 'just now';
    if (diff < 60) return `${diff}s ago`;
    if (diff < 3600) return `${Math.floor(diff / 60)}m ago`;
    return `${Math.floor(diff / 3600)}h ago`;
  });

  let overallPct = $derived(
    status?.pipeline_total > 0
      ? Math.round((status.pipeline_completed / status.pipeline_total) * 100)
      : 0
  );

  let isComplete = $derived(
    status && !status.pipeline_active && status.pipeline_running === 0
  );

  // Has any timing data?
  let hasTimeline = $derived(
    status?.timing && Object.keys(status.timing).length > 0
  );

  // Processes grouped by status
  let runningProcs = $derived(
    status ? Object.entries(status.processes)
      .filter(([, p]) => p.status === 'running')
      .sort(([a], [b]) => a.localeCompare(b)) : []
  );

  let errorProcs = $derived(
    status ? Object.entries(status.processes)
      .filter(([, p]) => p.status === 'failed' || p.status === 'warning')
      .sort(([, a], [, b]) => {
        // failed before warning
        if (a.status !== b.status) return a.status === 'failed' ? -1 : 1;
        return 0;
      }) : []
  );

  // --- Timeline / Gantt chart ---

  const LABEL_W = 160;
  const BAR_H = 16;
  const ROW_GAP = 4;
  const PAD_X = 16;
  const AXIS_H = 24;

  // Sort processes for timeline: by start time, then alphabetically
  let timelineProcs = $derived.by(() => {
    if (!status?.timing) return [];
    return Object.entries(status.timing)
      .map(([name, t]) => ({
        name,
        start: t.start,
        end: t.end,
        exit_code: t.exit_code,
        procStatus: status.processes[name]?.status,
      }))
      .sort((a, b) => {
        if (a.start && b.start) return a.start.localeCompare(b.start);
        if (a.start) return -1;
        if (b.start) return 1;
        return a.name.localeCompare(b.name);
      });
  });

  // Time bounds
  let tMin = $derived(status?.pipeline_start ? new Date(status.pipeline_start).getTime() : 0);
  let tMax = $derived.by(() => {
    if (status?.pipeline_active) return now;
    if (status?.pipeline_end) return new Date(status.pipeline_end).getTime();
    return now;
  });
  let tRange = $derived(Math.max(tMax - tMin, 1));

  // Chart dimensions
  let chartW = $state(500);
  let barAreaW = $derived(Math.max(chartW - LABEL_W - PAD_X * 2, 100));
  let svgH = $derived(AXIS_H + timelineProcs.length * (BAR_H + ROW_GAP) + 8);

  // Responsive width via container ref
  let containerEl = $state(null);
  $effect(() => {
    if (!containerEl) return;
    const ro = new ResizeObserver(entries => {
      for (const entry of entries) {
        chartW = entry.contentRect.width;
      }
    });
    ro.observe(containerEl);
    return () => ro.disconnect();
  });

  // Time axis ticks
  let ticks = $derived.by(() => {
    if (tRange <= 1) return [];
    const hours = tRange / 3600000;
    let intervalMs;
    if (hours < 1) intervalMs = 10 * 60000;       // 10 min
    else if (hours < 6) intervalMs = 30 * 60000;   // 30 min
    else if (hours < 24) intervalMs = 60 * 60000;  // 1 hour
    else intervalMs = 120 * 60000;                  // 2 hours

    const result = [];
    const firstTick = Math.ceil(tMin / intervalMs) * intervalMs;
    for (let t = firstTick; t <= tMax; t += intervalMs) {
      const x = LABEL_W + ((t - tMin) / tRange) * barAreaW;
      const d = new Date(t);
      const label = `${String(d.getHours()).padStart(2, '0')}:${String(d.getMinutes()).padStart(2, '0')}`;
      result.push({ x, label });
    }
    return result;
  });

  // Map time to X position
  function timeToX(isoStr) {
    if (!isoStr) return LABEL_W;
    const t = new Date(isoStr).getTime();
    return LABEL_W + ((t - tMin) / tRange) * barAreaW;
  }

  // Bar color based on process status / exit code
  function barColor(proc) {
    if (proc.procStatus === 'running') return '#f59e0b';   // amber
    if (proc.procStatus === 'failed' || proc.exit_code > 0) return '#f43f5e'; // rose
    if (proc.procStatus === 'warning') return '#f59e0b';   // amber
    return '#06b6d4';  // cyan
  }

  // Duration formatting for tooltip
  function formatDuration(startStr, endStr) {
    if (!startStr || !endStr) return '';
    const ms = new Date(endStr).getTime() - new Date(startStr).getTime();
    const sec = Math.floor(ms / 1000);
    if (sec < 60) return `${sec}s`;
    const min = Math.floor(sec / 60);
    if (min < 60) return `${min}m ${sec % 60}s`;
    const hr = Math.floor(min / 60);
    return `${hr}h ${min % 60}m`;
  }

  // Tooltip state
  let tooltip = $state(null);

  function showTooltip(proc, event) {
    const dur = proc.end ? formatDuration(proc.start, proc.end) : 'running...';
    const exitStr = proc.exit_code != null ? ` (exit ${proc.exit_code})` : '';
    const startTime = proc.start ? new Date(proc.start).toLocaleTimeString() : '';
    const endTime = proc.end ? new Date(proc.end).toLocaleTimeString() : '';
    tooltip = {
      text: `${proc.name}: ${dur}${exitStr}`,
      detail: proc.end ? `${startTime} \u2192 ${endTime}` : `Started ${startTime}`,
      x: event.clientX,
      y: event.clientY,
    };
  }

  function hideTooltip() {
    tooltip = null;
  }

  function formatEta(min) {
    if (min == null) return '';
    if (min < 60) return `~${Math.round(min)}m`;
    const h = min / 60;
    return `~${h.toFixed(1)}h`;
  }

  function formatElapsed(min) {
    if (min == null) return '';
    if (min < 60) return `${Math.round(min)}m`;
    const h = Math.floor(min / 60);
    const m = Math.round(min % 60);
    return m > 0 ? `${h}h ${m}m` : `${h}h`;
  }

  let copiedPath = $state(null);
  async function copyToClipboard(text) {
    try {
      await navigator.clipboard.writeText(text);
      copiedPath = text;
      setTimeout(() => { copiedPath = null; }, 2000);
    } catch {
      // Fallback: select text
    }
  }
</script>

{#if status}
  <div class="bg-slate-800 rounded-lg border border-slate-700 overflow-hidden mb-4">
    <!-- Header: overall progress (clickable when timeline data exists) -->
    <div
      class="px-4 py-3 border-b border-slate-700/50 {hasTimeline ? 'cursor-pointer select-none' : ''}"
      onclick={() => { if (hasTimeline) showTimeline = !showTimeline; }}
      onkeydown={(e) => { if (hasTimeline && (e.key === 'Enter' || e.key === ' ')) { e.preventDefault(); showTimeline = !showTimeline; }}}
      role={hasTimeline ? 'button' : undefined}
      tabindex={hasTimeline ? 0 : undefined}
    >
      <div class="flex items-center justify-between mb-2">
        <div class="flex items-center gap-3">
          {#if hasTimeline}
            <span class="text-xs text-slate-500 transition-transform duration-200 {showTimeline ? 'rotate-90' : ''}" style="display:inline-block">&#9654;</span>
          {/if}
          <h3 class="text-sm font-medium text-slate-300">
            {#if isComplete}
              Pipeline Complete
            {:else}
              Pipeline Progress
            {/if}
          </h3>
          <span class="text-xs text-slate-500">
            {status.pipeline_completed}/{status.pipeline_total} completed
          </span>
        </div>
        <div class="flex items-center gap-3 text-xs text-slate-500">
          {#if status.pipeline_failed > 0}
            <span class="text-rose-400">{status.pipeline_failed} failed</span>
          {/if}
          {#if status.pipeline_warning > 0}
            <span class="text-amber-400">{status.pipeline_warning} warnings</span>
          {/if}
          {#if updatedAgo}
            <span>Updated {updatedAgo}</span>
          {/if}
        </div>
      </div>

      <!-- Overall progress bar -->
      <div class="w-full bg-slate-700 rounded-full h-2 overflow-hidden">
        <div
          class="h-full rounded-full transition-all duration-700 ease-out {isComplete ? 'bg-cyan-500' : 'bg-cyan-500/80'}"
          style="width: {overallPct}%"
        ></div>
      </div>
      <div class="text-right text-xs text-slate-500 mt-1">{overallPct}%</div>
    </div>

    <!-- Timeline (Gantt chart) -->
    {#if showTimeline && hasTimeline}
      <div class="px-4 py-3 border-b border-slate-700/50 overflow-x-auto" bind:this={containerEl}>
        <svg width="100%" height={svgH} class="gantt-chart" style="min-width: 400px">
          <!-- Time axis -->
          <line x1={LABEL_W} y1={AXIS_H - 4} x2={LABEL_W + barAreaW} y2={AXIS_H - 4}
                stroke="#475569" stroke-width="1" />
          {#each ticks as tick}
            <line x1={tick.x} y1={AXIS_H - 8} x2={tick.x} y2={AXIS_H - 2}
                  stroke="#64748b" stroke-width="1" />
            <text x={tick.x} y={AXIS_H - 12} text-anchor="middle"
                  fill="#94a3b8" font-size="10" font-family="monospace">{tick.label}</text>
            <!-- Subtle gridline -->
            <line x1={tick.x} y1={AXIS_H} x2={tick.x} y2={svgH}
                  stroke="#334155" stroke-width="0.5" stroke-dasharray="2,4" />
          {/each}

          <!-- Process rows -->
          {#each timelineProcs as proc, i}
            {@const y = AXIS_H + i * (BAR_H + ROW_GAP)}
            {@const x1 = timeToX(proc.start)}
            {@const x2 = proc.end ? timeToX(proc.end) : timeToX(new Date(now).toISOString())}
            {@const barW = Math.max(x2 - x1, 3)}
            {@const color = barColor(proc)}

            <!-- Process name label -->
            <text x={LABEL_W - 8} y={y + BAR_H / 2 + 4} text-anchor="end"
                  fill="#94a3b8" font-size="11" font-family="monospace">{proc.name}</text>

            <!-- Bar -->
            <rect
              x={x1} y={y} width={barW} height={BAR_H} rx="2"
              fill={color} opacity="0.8"
              onmouseenter={(e) => showTooltip(proc, e)}
              onmouseleave={hideTooltip}
            >
              {#if proc.procStatus === 'running'}
                <animate attributeName="opacity" values="0.8;0.4;0.8" dur="2s" repeatCount="indefinite" />
              {/if}
            </rect>
          {/each}
        </svg>
      </div>
    {/if}

    <!-- Running processes -->
    {#if runningProcs.length > 0}
      <div class="px-4 py-3 border-b border-slate-700/50">
        {#each runningProcs as [name, proc]}
          <div class="mb-3 last:mb-0">
            <div class="flex items-center justify-between mb-1">
              <div class="flex items-center gap-2">
                <span class="inline-block w-2 h-2 rounded-full bg-amber-400 animate-pulse"></span>
                <span class="text-xs font-medium text-amber-300">{name}</span>
              </div>
              <div class="flex items-center gap-3 text-xs text-slate-500">
                {#if proc.progress != null}
                  <span class="text-amber-300/80">{Math.round(proc.progress * 100)}%</span>
                {/if}
                {#if proc.eta_min != null}
                  <span>ETA {formatEta(proc.eta_min)}</span>
                {/if}
                {#if proc.elapsed_min != null && proc.progress == null}
                  <span>running for {formatElapsed(proc.elapsed_min)}</span>
                {/if}
              </div>
            </div>

            {#if proc.progress != null}
              <div class="w-full bg-slate-700 rounded-full h-1.5 overflow-hidden mb-1">
                <div
                  class="h-full rounded-full bg-amber-400/60 transition-all duration-700"
                  style="width: {Math.round(proc.progress * 100)}%"
                ></div>
              </div>
            {/if}

            {#if proc.progress_detail}
              <p class="text-xs text-slate-500 ml-4">{proc.progress_detail}</p>
            {/if}
          </div>
        {/each}
      </div>
    {/if}

    <!-- Errors and warnings -->
    {#if errorProcs.length > 0}
      <div class="px-4 py-3">
        {#each errorProcs as [name, proc]}
          <div class="mb-3 last:mb-0">
            <div class="flex items-start gap-2">
              <span class="mt-0.5 text-sm {proc.status === 'failed' ? 'text-rose-400' : 'text-amber-400'}">
                {proc.status === 'failed' ? '\u2716' : '\u26A0'}
              </span>
              <div class="flex-1 min-w-0">
                <div class="flex items-center gap-2 mb-0.5">
                  <span class="text-xs font-medium {proc.status === 'failed' ? 'text-rose-300' : 'text-amber-300'}">{name}</span>
                  {#if proc.exit_code != null}
                    <span class="text-xs text-slate-600">exit {proc.exit_code}</span>
                  {/if}
                </div>

                {#if proc.error}
                  <p class="text-xs text-slate-400 break-words">{proc.error}</p>
                {/if}

                {#if proc.work_dir}
                  <div class="mt-1.5">
                    <button
                      onclick={() => copyToClipboard(proc.work_dir)}
                      class="text-xs font-mono text-slate-600 hover:text-slate-400 transition-colors cursor-pointer break-all text-left"
                      title="Copy work directory path"
                    >
                      {proc.work_dir}{copiedPath === proc.work_dir ? ' \u2713' : ''}
                    </button>
                    {#if proc.files}
                      <div class="flex gap-2 mt-1">
                        {#each Object.entries(proc.files) as [label, filename]}
                          <button
                            onclick={() => copyToClipboard(proc.work_dir + '/' + filename)}
                            class="text-xs font-mono text-slate-600 hover:text-slate-400 transition-colors cursor-pointer"
                            title="Copy path: {filename}"
                          >
                            {filename}{copiedPath === proc.work_dir + '/' + filename ? ' \u2713' : ''}
                          </button>
                        {/each}
                      </div>
                    {/if}
                  </div>
                {/if}
              </div>
            </div>
          </div>
        {/each}
      </div>
    {/if}
  </div>
{/if}

<!-- Floating tooltip -->
{#if tooltip}
  <div
    class="fixed z-50 px-2 py-1 bg-slate-900 border border-slate-600 rounded text-xs text-slate-200 pointer-events-none shadow-lg"
    style="left: {tooltip.x + 12}px; top: {tooltip.y - 8}px"
  >
    <div class="font-medium">{tooltip.text}</div>
    <div class="text-slate-400">{tooltip.detail}</div>
  </div>
{/if}

<style>
  @keyframes pulse {
    0%, 100% { opacity: 1; }
    50% { opacity: 0.4; }
  }
  .animate-pulse {
    animation: pulse 2s ease-in-out infinite;
  }
  .gantt-chart rect {
    cursor: pointer;
  }
  .gantt-chart rect:hover {
    opacity: 1 !important;
    filter: brightness(1.2);
  }
</style>
