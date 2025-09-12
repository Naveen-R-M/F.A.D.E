'use client';

import TopBar from '@/components/header/TopBar';
import SegmentedTabs from '@/components/nav/SegmentedTabs';
import ChatUI from '@/components/chat/ChatUI';
import InProgressPanel from '@/components/progress/InProgressPanel';
import HistoryPanel from '@/components/history/HistoryPanel';
import { useEffect, useState } from 'react';
import { ModelProvider } from '@/components/context/ModelContext';
import { JobsProvider } from '@/components/context/JobsContext';
import { ChatProvider } from '@/components/context/ChatContext';

type Tab = 'chat' | 'history' | 'progress';

export default function Page() {
  const [tab, setTab] = useState<Tab>('chat');

  // (optional) remember last tab in session
  useEffect(() => {
    const t = sessionStorage.getItem('fade_active_tab') as Tab | null;
    if (t) setTab(t);
  }, []);
  useEffect(() => {
    sessionStorage.setItem('fade_active_tab', tab);
  }, [tab]);

  return (
    <ModelProvider>
      <JobsProvider>
        <ChatProvider>
          <main className="mx-auto max-w-4xl px-4 py-4">
            <TopBar />
            <div className="mt-4 flex justify-center">
              <SegmentedTabs value={tab} onChange={setTab} size="sm" />
            </div>

            <section className="mt-10">
              {tab === 'chat' && <ChatUI />}
              {tab === 'history' && <HistoryPanel onOpen={() => setTab('chat')} />}
              {tab === 'progress' && <InProgressPanel />}
            </section>
          </main>
        </ChatProvider>
      </JobsProvider>
    </ModelProvider>
  );
}
