'use client';
import { useModel, ModelId } from '@/components/context/ModelContext';

const OPTIONS: { id: ModelId; label: string }[] = [
  { id: 'gemini-2.5-pro',  label: 'Gemini 2.5 Pro' },
  { id: 'gemini-2.5-flash',  label: 'Gemini 2.5 Flash' },
  { id: 'gemini-2.5-flash-lite',  label: 'Gemini 2.5 Flash Lite' },
  { id: 'local',   label: 'Local' },
];

export default function ModelSelector() {
  const { model, setModel } = useModel();

  return (
    <label className="inline-flex items-center gap-2">
      <span className="sr-only">Select model</span>
      <select
        value={model}
        onChange={(e) => setModel(e.target.value as ModelId)}
        className="rounded-lg border border-white/10 bg-white/5 px-3 py-2 text-sm text-white outline-none hover:bg-white/10"
      >
        {OPTIONS.map((m) => (
          <option key={m.id} value={m.id} className="bg-[#0b0b0b]">
            {m.label}
          </option>
        ))}
      </select>
    </label>
  );
}
