export default function AnimatedOrb({ gradient }: { gradient: string }) {
  return (
    <div className="relative h-14 w-14 animate-float">
      {/* soft glow */}
      <div className={`absolute inset-0 rounded-full ${gradient} blur-xl opacity-70 animate-pulse-glow`} />
      {/* core */}
      <div className={`absolute inset-0 rounded-full ${gradient}`} />
      {/* rotating highlight ring */}
      <div className="absolute inset-0 rounded-full p-[2px]">
        <div className="h-full w-full rounded-full bg-[conic-gradient(from_0deg,rgba(255,255,255,.35),transparent_60%,transparent_100%)] animate-spin-slow" />
      </div>
      {/* subtle inner shine */}
      <div className="absolute inset-0 rounded-full bg-[radial-gradient(circle_at_35%_30%,rgba(255,255,255,.35),transparent_55%)] mix-blend-screen" />
    </div>
  );
}
