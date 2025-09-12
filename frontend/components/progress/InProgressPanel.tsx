// components/progress/InProgressPanel.tsx
'use client';
import { useEffect, useMemo, useState } from 'react';
import { useJobs } from '@/components/context/JobsContext';

type Filter = 'all' | 'queued' | 'running' | 'completed' | 'failed';

export default function InProgressPanel() {
  const { jobs, deleteJob } = useJobs(); // ✅ pull delete
  const [filter, setFilter] = useState<Filter>('all');

  useEffect(() => {
    const saved = sessionStorage.getItem('fade_jobs_filter') as Filter | null;
    if (saved) setFilter(saved);
  }, []);
  useEffect(() => {
    sessionStorage.setItem('fade_jobs_filter', filter);
  }, [filter]);

  const counts = useMemo(() => {
    const base = { all: jobs.length, queued: 0, running: 0, completed: 0, failed: 0 } as Record<Filter, number>;
    for (const j of jobs) base[j.status as Filter] = (base[j.status as Filter] ?? 0) + 1;
    return base;
  }, [jobs]);

  const visible = useMemo(() => (filter === 'all' ? jobs : jobs.filter(j => j.status === filter)), [jobs, filter]);

  return (
    <div className="mx-auto max-w-3xl space-y-3">
      <StatusFilter filter={filter} setFilter={setFilter} counts={counts} />

      {visible.length === 0 ? (
        <div className="rounded-2xl border border-white/10 bg-white/5 p-6 text-center text-sm text-white/75">
          {filter === 'all' ? 'No jobs yet.' : `No ${filter} jobs.`}
        </div>
      ) : (
        visible.map((j) => (
          <div key={j.id} className="rounded-xl border border-white/10 bg-white/5 p-4 text-sm">
            <div className="flex items-start justify-between gap-3">
              <div className="min-w-0 flex-1">
                <div className="flex items-center justify-between gap-3">
                  <div className="font-medium truncate">{j.query}</div>
                  <StatusBadge status={j.status} />
                </div>
                <div className="mt-1 text-white/60">
                  <span className="mr-3">Model: {j.model}</span>
                  <span className="mr-3">User: {j.user_id.slice(0, 8)}…</span>
                  <span>{new Date(j.current_time).toLocaleString()}</span>
                </div>
                {j.error && <div className="mt-2 text-rose-300">{j.error}</div>}
              </div>

              {/* ✅ Delete button (same style vibe as History) */}
              <button
                type="button"
                aria-label="Delete job"
                title="Delete"
                onClick={(e) => {
                  e.stopPropagation();
                  const ok = window.confirm('Remove this job from the list?');
                  if (ok) deleteJob(j.id);
                }}
                className="shrink-0 rounded-full p-2 text-white/80 hover:bg-white/10 hover:text-white"
              >
                <svg viewBox="0 0 24 24" className="h-4 w-4" fill="none" stroke="currentColor" strokeWidth="1.6">
                  <path d="M4 7h16" />
                  <path d="M10 11v6M14 11v6" />
                  <path d="M6 7l1 13a2 2 0 0 0 2 2h6a2 2 0 0 0 2-2l1-13" />
                  <path d="M9 7V5a2 2 0 0 1 2-2h2a2 2 0 0 1 2 2v2" />
                </svg>
              </button>
            </div>
          </div>
        ))
      )}
    </div>
  );
}

function StatusFilter({
  filter,
  setFilter,
  counts,
}: {
  filter: Filter;
  setFilter: (f: Filter) => void;
  counts: Record<Filter, number>;
}) {
  const items: { key: Filter; label: string }[] = [
    { key: 'all', label: 'All' },
    { key: 'queued', label: 'Queued' },
    { key: 'running', label: 'Running' },
    { key: 'completed', label: 'Completed' },
    { key: 'failed', label: 'Failed' },
  ];

  return (
    <div
      role="tablist"
      aria-label="Job status filter"
      className="inline-flex items-center rounded-full border border-white/10 bg-white/5 p-1"
    >
      {items.map((it) => {
        const active = it.key === filter;
        return (
          <button
            key={it.key}
            type="button"
            role="tab"
            aria-selected={active}
            onClick={() => setFilter(it.key)}
            className={[
              'relative rounded-full px-3 py-1.5 text-xs transition-colors',
              active ? 'bg-indigo-500/30 text-white' : 'text-white/85 hover:bg-white/7',
            ].join(' ')}
          >
            <span className="align-middle">{it.label}</span>
            <span className="ml-2 rounded-full bg-white/10 px-1.5 py-0.5 text-[10px]">
              {counts[it.key]}
            </span>
          </button>
        );
      })}
    </div>
  );
}

function StatusBadge({
  status,
}: {
  status: 'queued' | 'running' | 'completed' | 'failed';
}) {
  const style =
    status === 'queued'
      ? 'bg-amber-500/20 text-amber-300'
      : status === 'running'
      ? 'bg-sky-500/20 text-sky-300'
      : status === 'completed'
      ? 'bg-emerald-500/20 text-emerald-300'
      : 'bg-rose-500/20 text-rose-300';
  return <span className={`rounded-full px-2 py-0.5 text-xs capitalize ${style}`}>{status}</span>;
}
