'use client';
import { createContext, useContext, useState } from 'react';

export type ModelId = 'gpt-4o' | 'gpt-4.1' | 'gpt-3.5' | 'local';

type Ctx = {
  model: ModelId;
  setModel: (m: ModelId) => void;
};

const ModelContext = createContext<Ctx | null>(null);

export function ModelProvider({ children }: { children: React.ReactNode }) {
  const [model, setModel] = useState<ModelId>('gpt-4o');
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
