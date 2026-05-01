<script>
  // Feedback drawer — ports microscape-app's FeedbackForm.svelte UI/payload.
  // Always submits to https://microscape.app/api/feedback so feedback from
  // standalone viz deployments lands in the same backend queue.

  let { open = false, onClose = () => {} } = $props();

  const ENDPOINT = 'https://microscape.app/api/feedback';

  let message = $state('');
  let sending = $state(false);
  let sent = $state(false);
  let error = $state('');

  async function submit() {
    if (!message.trim()) return;
    sending = true;
    error = '';
    try {
      const pageUrl = (typeof window !== 'undefined')
        ? window.location.href
        : '';
      const res = await fetch(ENDPOINT, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ page_url: pageUrl, message: message.trim() })
      });
      if (res.ok) {
        sent = true;
        message = '';
        setTimeout(() => {
          sent = false;
          onClose();
        }, 1800);
      } else {
        let body = null;
        try { body = await res.json(); } catch { /* not JSON */ }
        error = (body && body.error) || `Server returned ${res.status}`;
      }
    } catch (e) {
      error = 'Could not reach feedback service. Please try again later.';
    } finally {
      sending = false;
    }
  }

  function handleClose() {
    message = '';
    error = '';
    sent = false;
    onClose();
  }
</script>

{#if open}
  <div class="fixed right-0 top-14 bottom-0 w-96 bg-slate-900 border-l border-slate-700 z-40 flex flex-col shadow-2xl">
    <div class="p-4 border-b border-slate-700 flex items-center justify-between">
      <h3 class="text-sm font-semibold text-slate-200">Feedback</h3>
      <button
        class="text-xs px-2 py-1 rounded text-slate-400 hover:text-slate-200 border border-slate-600 hover:border-slate-500 transition-colors"
        onclick={handleClose}
      >Close</button>
    </div>
    <div class="flex-1 overflow-y-auto p-4 space-y-3">
      <p class="text-xs text-slate-400">
        Report an issue or suggestion. Submissions are sent to the microscape.app team.
      </p>
      {#if sent}
        <p class="text-sm text-emerald-400">Thanks! Your feedback has been recorded.</p>
      {:else}
        <form onsubmit={(e) => { e.preventDefault(); submit(); }} class="space-y-3">
          <textarea
            bind:value={message}
            placeholder="What's wrong or what could be better?"
            rows="6"
            class="w-full px-3 py-2 bg-slate-800 border border-slate-700 rounded-md text-slate-100 text-sm placeholder-slate-500 focus:outline-none focus:border-cyan-500 resize-y"
          ></textarea>
          {#if error}
            <p class="text-xs text-rose-400">{error}</p>
          {/if}
          <div class="flex justify-end gap-2">
            <button
              type="button"
              onclick={handleClose}
              class="px-3 py-1.5 text-sm text-slate-400 hover:text-slate-200 transition-colors"
            >Cancel</button>
            <button
              type="submit"
              disabled={sending || !message.trim()}
              class="px-3 py-1.5 bg-cyan-600 text-white rounded-md hover:bg-cyan-500 disabled:opacity-50 disabled:cursor-not-allowed text-sm font-medium transition-colors"
            >{sending ? 'Sending…' : 'Send feedback'}</button>
          </div>
        </form>
      {/if}
    </div>
  </div>
{/if}
