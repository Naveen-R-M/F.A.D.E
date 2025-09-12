export function UserAvatar({ size = 28 }: { size?: number }) {
  return (
    <div
      style={{ width: size, height: size }}
      className="grid place-items-center rounded-full bg-gradient-to-br from-blue-500 to-cyan-400 text-[10px] font-semibold text-black"
    >
      {/* Replace with user image later */}
      <span className="sr-only">User profile</span>
      <span className="pointer-events-none select-none">•</span>
    </div>
  );
}