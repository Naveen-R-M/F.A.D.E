export function BellIcon({ className = '' }: { className?: string }) {
  return (
    <svg viewBox="0 0 24 24" fill="none" className={className} aria-hidden="true">
      <path
        d="M6 9a6 6 0 1112 0c0 4 2 5 2 5H4s2-1 2-5z"
        stroke="currentColor"
        strokeWidth="1.5"
      />
      <path d="M9 18a3 3 0 006 0" stroke="currentColor" strokeWidth="1.5" />
    </svg>
  );
}
