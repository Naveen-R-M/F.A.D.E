'use client';

import { useId, useState } from 'react';

type TabKey = 'chat' | 'history' | 'progress';

const TABS: { key: TabKey; label: string }[] = [
  { key: 'chat',    label: 'Chat' },
  { key: 'history', label: 'History' },
  { key: 'progress', label: 'In Progress' },
];

export default function SegmentedTabs({
  value,
  onChange,
  className = '',
  size = 'sm',
}: {
  value?: TabKey;
  onChange?: (v: TabKey) => void;
  className?: string;
  /** 'sm' | 'md' */
  size?: 'sm' | 'md';
}) {
  const [internal, setInternal] = useState<TabKey>(value ?? 'chat');
  const selected = value ?? internal;
  const setSelected = (v: TabKey) => {
    setInternal(v);
    onChange?.(v);
  };

  const group = useId();

  const sizes =
    size === 'md'
      ? { container: 'p-1.5', tab: 'px-4 py-2 text-sm', dot: 'h-2 w-2' }
      : { container: 'p-1', tab: 'px-3 py-1.5 text-xs', dot: 'h-1.5 w-1.5' };

  return (
    <div
      role="tablist"
      aria-label="Primary"
      className={`inline-flex items-center rounded-full border border-white/10 bg-white/5 shadow-sm ${sizes.container} ${className}`}
      onKeyDown={(e) => {
        const idx = TABS.findIndex(t => t.key === selected);
        if (e.key === 'ArrowRight') {
          e.preventDefault();
          setSelected(TABS[(idx + 1) % TABS.length].key);
        } else if (e.key === 'ArrowLeft') {
          e.preventDefault();
          setSelected(TABS[(idx - 1 + TABS.length) % TABS.length].key);
        }
      }}
    >
      {TABS.map((t) => {
        const active = t.key === selected;
        return (
          <button
            key={t.key}
            id={`${group}-${t.key}`}
            role="tab"
            aria-selected={active}
            aria-controls={`${group}-${t.key}-panel`}
            className={[
              'relative inline-flex items-center gap-2 rounded-full transition-colors',
              sizes.tab,
              active
                ? 'bg-indigo-500/30 text-white shadow-[inset_0_0_0_1px_rgba(255,255,255,0.1)]'
                : 'text-white/85 hover:bg-white/7'
            ].join(' ')}
            onClick={() => setSelected(t.key)}
          >
            <span className="font-medium">{t.label}</span>
          </button>
        );
      })}
    </div>
  );
}
