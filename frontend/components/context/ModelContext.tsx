'use client';
import { createContext, useContext, useState } from 'react';

// Updated to match Gemini models from constants
export type ModelId = 'gemini-2.5-pro' | 'gemini-2.5-flash' | 'gemini-2.5-flash-lite';

type Ctx = {
  model: ModelId;
  setModel: (m: ModelId) => void;
};

const ModelContext = createContext<Ctx | null>(null);

export function ModelProvider({ children }: { children: React.ReactNode }) {
  const [model, setModel] = useState<ModelId>('gemini-2.5-pro');
  return (
    <ModelContext.Provider value={{ model, setModel }}>
      {children}
    </ModelContext.Provider>
  );
}

export function useModel() {
  const ctx = useContext(ModelContext);
  if (!ctx) throw new Error('useModel must be used within <ModelProvider>');
  return ctx;
}