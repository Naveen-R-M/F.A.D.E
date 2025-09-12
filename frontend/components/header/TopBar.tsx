'use client';

import ModelSelector from './ModelSelector';
import { BellIcon } from '@/components/icons/BellIcon';
import { UserAvatar } from '@/components/icons/UserAvatar';
import Logo from './Logo';

export default function TopBar() {
  return (
    <header
      className={
        // glassmorphism: translucent background + backdrop blur, thin border and subtle shadow
        "relative w-full sticky top-0 z-40 rounded-xl overflow-hidden px-4 md:px-7 py-3 bg-white/5 backdrop-blur-md backdrop-saturate-150 border-b border-white/10 shadow-sm"
      }
    >
      {/* blurred border layer for smoothness */}
      <div
        aria-hidden
        className="absolute inset-0 -z-10 rounded-xl bg-white/5 blur-sm pointer-events-none"
      />
      {/* 3 equal columns => true center logo regardless of side widths */}
      <div
        className="
          grid grid-cols-3 items-center
        "
      >
        {/* Left: model selector flush toward left margin */}
        <div className="justify-self-start">
          <ModelSelector />
        </div>

        {/* Center: logo perfectly centered */}
        <div className="justify-self-center">
          <Logo text="F.A.D.E" sizeRem={1.4} />
        </div>

        {/* Right: actions flush to right margin */}
        <nav className="justify-self-end flex items-center gap-2">
          <IconButton ariaLabel="Notifications" badge>
            <BellIcon className="h-5 w-5" />
          </IconButton>
          <button
            aria-label="Profile"
            className="ml-1 inline-flex items-center gap-2 rounded-full p-1 hover:bg-white/5"
          >
            <UserAvatar size={30} />
          </button>
        </nav>
      </div>
    </header>
  );
}

function IconButton({
  children,
  ariaLabel,
  badge = false,
}: {
  children: React.ReactNode;
  ariaLabel: string;
  badge?: boolean;
}) {
  return (
    <button
      aria-label={ariaLabel}
      className="relative inline-flex items-center justify-center rounded-full p-2 hover:bg-white/5 focus:outline-none focus-visible:ring-2 focus-visible:ring-white/30"
    >
      {children}
      {badge && (
        <span className="absolute right-1 top-1 block h-2 w-2 rounded-full bg-orange-500" />
      )}
    </button>
  );
}
