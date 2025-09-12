// components/context/JobsContext.tsx
'use client';

import { createContext, useContext, useEffect, useMemo, useRef, useState } from 'react';

export type JobStatus = 'queued' | 'running' | 'completed' | 'failed';
export type Job = {
  id: string;
  query: string;
  model: string;
  user_id: string;
  current_time: string; // ISO
  status: JobStatus;
  error?: string;
};

type Ctx = {
  jobs: Job[];
  addJob: (j: Job) => void;
  updateJob: (id: string, patch: Partial<Job>) => void;
  clearCompleted: () => void;
  clearAll: () => void;
  deleteJob: (id: string) => void;
};
const JobsContext = createContext<Ctx | null>(null);
const STORAGE_KEY = 'fade_jobs_v1';

const isBrowser = () => typeof window !== 'undefined';

export function JobsProvider({ children }: { children: React.ReactNode }) {
  const [jobs, setJobs] = useState<Job[]>([]);

  // Load
  useEffect(() => {
    if (!isBrowser()) return;
    try {
      const raw = localStorage.getItem(STORAGE_KEY);
      if (raw) {
        const parsed = JSON.parse(raw) as Job[];
        if (Array.isArray(parsed)) setJobs(parsed);
      }
    } catch { }
  }, []);

  // Persist (debounced for normal updates)
  const saveTimer = useRef<number | null>(null);
  const save = (payload: Job[]) => {
    try {
      localStorage.setItem(STORAGE_KEY, JSON.stringify(payload));
      // cross-tab notify
      window.dispatchEvent(new StorageEvent('storage', { key: STORAGE_KEY, newValue: JSON.stringify(payload) }));
    } catch { }
  };
  useEffect(() => {
    if (!isBrowser()) return;
    if (saveTimer.current) window.clearTimeout(saveTimer.current);
    saveTimer.current = window.setTimeout(() => save(jobs), 120);
    return () => {
      if (saveTimer.current) window.clearTimeout(saveTimer.current);
    };
  }, [jobs]);

  // Cross-tab sync
  useEffect(() => {
    if (!isBrowser()) return;
    const onStorage = (e: StorageEvent) => {
      if (e.key === STORAGE_KEY && e.newValue) {
        try {
          const parsed = JSON.parse(e.newValue) as Job[];
          if (Array.isArray(parsed)) setJobs(parsed);
        } catch { }
      }
    };
    window.addEventListener('storage', onStorage);
    return () => window.removeEventListener('storage', onStorage);
  }, []);

  // API
  const addJob = (j: Job) => setJobs(prev => [{ ...j, status: j.status.toLowerCase() as JobStatus }, ...prev]);

  const updateJob = (id: string, patch: Partial<Job>) =>
    setJobs(prev => prev.map(job => job.id === id ? { ...job, ...patch, status: (patch.status ? String(patch.status).toLowerCase() : job.status) as JobStatus } : job));

  const clearCompleted = () => setJobs(prev => prev.filter(j => j.status !== 'completed'));

  const clearAll = () => setJobs([]);

  const deleteJob = (id: string) => setJobs(prev => prev.filter(j => j.id !== id)); // ✅ remove one

  const value: Ctx = { jobs, addJob, updateJob, clearCompleted, clearAll, deleteJob }; // ✅ expose
  return <JobsContext.Provider value={value}>{children}</JobsContext.Provider>;
}

export function useJobs() {
  const ctx = useContext(JobsContext);
  if (!ctx) throw new Error('useJobs must be used within <JobsProvider>');
  return ctx;
}
