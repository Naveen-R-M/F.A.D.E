// components/context/ChatContext.tsx - Enhanced to support message details
'use client';

import React, {
  createContext,
  useContext,
  useEffect,
  useMemo,
  useRef,
  useState,
  type ReactNode,
} from 'react';

export type Role = 'user' | 'assistant' | 'system';

// Enhanced message type to support detailed logs
export type ChatMessage = { 
  id: string; 
  role: Role; 
  content: string; 
  ts: string;
  details?: any[]; // Store SSE events for "Show Your Work"
  isThinking?: boolean; // For thinking state
  currentStep?: string; // Current processing step
};

export type Conversation = {
  id: string;
  title: string;
  model: string;
  createdAt: string;
  updatedAt: string;
  messages: ChatMessage[];
};

type Ctx = {
  conversations: Conversation[];
  activeId: string | null;
  active?: Conversation;
  createConversation: (opts?: { model?: string; greeting?: string }) => string;
  setActive: (id: string) => void;
  appendMessage: (msg: Omit<ChatMessage, 'id' | 'ts'>, opts?: { convoId?: string }) => void;
  updateTitle: (id: string, title: string) => void;
  clearAll: () => void;
  deleteConversation: (id: string) => void;
};

const ChatContext = createContext<Ctx | null>(null);
const STORAGE_KEY = 'fade_conversations_v1';

export function ChatProvider({ children }: { children: ReactNode }) {
  const [conversations, setConversations] = useState<Conversation[]>([]);
  const [activeId, setActiveId] = useState<string | null>(null);

  // Load once
  useEffect(() => {
    try {
      const raw = localStorage.getItem(STORAGE_KEY);
      if (raw) {
        const parsed = JSON.parse(raw) as { conversations: Conversation[]; activeId: string | null };
        setConversations(parsed.conversations ?? []);
        setActiveId(parsed.activeId ?? null);
      }
    } catch {}
  }, []);

  // Persist on changes
  useEffect(() => {
    localStorage.setItem(STORAGE_KEY, JSON.stringify({ conversations, activeId }));
  }, [conversations, activeId]);

  // Ensure activeId is valid when list changes
  useEffect(() => {
    if (conversations.length === 0) {
      if (activeId !== null) setActiveId(null);
      return;
    }
    if (!conversations.some((c) => c.id === activeId)) {
      setActiveId(conversations[0].id);
    }
  }, [conversations, activeId]);

  const createConversation: Ctx['createConversation'] = (opts) => {
    const id = crypto.randomUUID();
    const now = new Date().toISOString();
    const model = opts?.model ?? 'gemini-2.5-pro';
    const greeting = opts?.greeting ?? 'Hi! Ask me anything.';
    const convo: Conversation = {
      id,
      title: 'New chat',
      model,
      createdAt: now,
      updatedAt: now,
      messages: [{ id: crypto.randomUUID(), role: 'assistant', content: greeting, ts: now }],
    };
    setConversations((prev) => [convo, ...prev]);
    setActiveId(id);
    return id;
  };

  const deleteConversation: Ctx['deleteConversation'] = (id) => {
    setConversations((prev) => {
      const next = prev.filter((c) => c.id !== id);
      if (id === activeId || !next.some((c) => c.id === activeId)) {
        setActiveId(next[0]?.id ?? null);
      }
      return next;
    });
  };

  const appendMessage: Ctx['appendMessage'] = (msg, opts) => {
    let targetId = opts?.convoId ?? activeId;

    // Lazily create if none is active
    if (!targetId) {
      const id = crypto.randomUUID();
      const now = new Date().toISOString();
      const convo: Conversation = {
        id,
        title: msg.role === 'user' ? msg.content.slice(0, 60) : 'New chat',
        model: 'gemini-2.5-pro',
        createdAt: now,
        updatedAt: now,
        messages: [],
      };
      setConversations((prev) => [convo, ...prev]);
      setActiveId(id);
      targetId = id;
    }

    const now = new Date().toISOString();
    setConversations((prev) =>
      prev.map((c) =>
        c.id !== targetId
          ? c
          : {
              ...c,
              updatedAt: now,
              title:
                c.title === 'New chat' && msg.role === 'user'
                  ? msg.content.slice(0, 60)
                  : c.title,
              messages: [...c.messages, { id: crypto.randomUUID(), ts: now, ...msg }],
            }
      )
    );
  };

  const updateTitle: Ctx['updateTitle'] = (id, title) =>
    setConversations((prev) => prev.map((c) => (c.id === id ? { ...c, title } : c)));

  const clearAll: Ctx['clearAll'] = () => {
    setConversations([]);
    setActiveId(null);
    localStorage.removeItem(STORAGE_KEY);
  };

  const value = useMemo<Ctx>(
    () => ({
      conversations,
      activeId,
      active: conversations.find((c) => c.id === activeId),
      createConversation,
      setActive: setActiveId,
      appendMessage,
      updateTitle,
      clearAll,
      deleteConversation,
    }),
    [conversations, activeId]
  );

  return <ChatContext.Provider value={value}>{children}</ChatContext.Provider>;
}

export function useChat() {
  const ctx = useContext(ChatContext);
  if (!ctx) throw new Error('useChat must be used within <ChatProvider>');
  return ctx;
}