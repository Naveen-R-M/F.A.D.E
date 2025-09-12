'use client';

import { useChat } from '@/components/context/ChatContext';

export default function HistoryPanel({ onOpen }: { onOpen: () => void }) {
  const { conversations, setActive, deleteConversation } = useChat();

  if (!conversations.length) {
    return (
      <div className="mx-auto max-w-3xl rounded-2xl border border-white/10 bg-white/5 p-6 text-center text-sm text-white/75">
        No conversations yet.
      </div>
    );
  }

  return (
    <div className="mx-auto max-w-3xl space-y-3">
      {conversations
        .slice()
        .sort((a, b) => +new Date(b.updatedAt) - +new Date(a.updatedAt))
        .map((c) => (
          <div
            key={c.id}
            className="flex items-start justify-between rounded-xl border border-white/10 bg-white/5 p-4"
          >
            {/* Clickable area to open the conversation */}
            <button
              type="button"
              onClick={() => { setActive(c.id); onOpen(); }}
              className="min-w-0 flex-1 text-left hover:bg-white/5 rounded-lg px-2 py-1 transition-colors"
              title="Open conversation"
            >
              <div className="flex items-center justify-between">
                <div className="truncate font-medium">{c.title || 'Untitled chat'}</div>
                <div className="shrink-0 text-xs text-white/60">
                  {new Date(c.updatedAt).toLocaleString()}
                </div>
              </div>
              <div className="mt-1 line-clamp-1 text-sm text-white/60">
                {c.messages.at(-1)?.content ?? ''}
              </div>
            </button>

            {/* Delete action */}
            <button
              type="button"
              aria-label="Delete conversation"
              className="ml-3 rounded-full p-2 text-white/80 hover:bg-white/10 hover:text-white"
              title="Delete"
              onClick={(e) => {
                e.stopPropagation(); // don't trigger open
                const ok = window.confirm('Delete this conversation? This cannot be undone.');
                if (ok) deleteConversation(c.id);
              }}
            >
              {/* tiny trash icon (SVG, no deps) */}
              <svg viewBox="0 0 24 24" className="h-4 w-4" fill="none" stroke="currentColor" strokeWidth="1.6">
                <path d="M4 7h16" />
                <path d="M10 11v6M14 11v6" />
                <path d="M6 7l1 13a2 2 0 0 0 2 2h6a2 2 0 0 0 2-2l1-13" />
                <path d="M9 7V5a2 2 0 0 1 2-2h2a2 2 0 0 1 2 2v2" />
              </svg>
            </button>
          </div>
        ))}
    </div>
  );
}
