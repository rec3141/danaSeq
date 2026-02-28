<script>
  let { status = null } = $props();

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
    <!-- Header: overall progress -->
    <div class="px-4 py-3 border-b border-slate-700/50">
      <div class="flex items-center justify-between mb-2">
        <div class="flex items-center gap-3">
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

<style>
  @keyframes pulse {
    0%, 100% { opacity: 1; }
    50% { opacity: 0.4; }
  }
  .animate-pulse {
    animation: pulse 2s ease-in-out infinite;
  }
</style>
