// components/header/Logo.tsx
'use client';

import { geistSans } from '@/lib/fonts'; // or import { Geist } from "next/font/google" and create here

type Props = {
  text?: string;        // default "F.A.D.E"
  sizeRem?: number;     // wordmark size
  className?: string;
};

export default function Logo({ text = 'F.A.D.E', sizeRem = 1.35, className = '' }: Props) {
  // We’ll do a Next.js-like treatment:
  //  - render chars individually to massage spacing
  //  - shrink/offset '.' dots
  //  - add a diagonal slash over the first 'A'
  const chars = text.toUpperCase().split('');

  return (
    <div
      className={`${geistSans.className} ${className} select-none leading-none`}
      style={{
        fontSize: `${sizeRem}rem`,
        fontWeight: 900,                 // Geist Black
        letterSpacing: '-0.02em',        // tight like Next.js
      }}
      aria-label={`${text} logo`}
    >
      {chars.map((ch, i) => {
        const isDot = ch === '.';
        const isA   = ch === 'A';

        if (isDot) {
          // small, lowered dot (Next.js vibe for ".js")
          return (
            <span
              key={`${ch}-${i}`}
              className="inline-block align-baseline"
              style={{
                transform: 'translateY(0.18em)',
                fontSize: '0.45em',
                marginLeft: '0.04em',
                marginRight: '0.04em',
                opacity: 0.9,
              }}
            >
              •
            </span>
          );
        }

        if (isA) {
          // diagonal slash through the A
          return (
            <span key={`${ch}-${i}`} className="relative inline-block">
              <span
                className="inline-block"
                style={{
                  paddingLeft: i === 0 ? '0.01em' : '0', // optical left edge
                  paddingRight: '0.01em',
                }}
              >
                A
              </span>
              <span
                aria-hidden
                className="pointer-events-none absolute inset-0"
                style={{
                  // a thin line cutting across the A, similar angle to Next.js
                  transform: 'rotate(70deg)',
                }}
              >
                <span
                  className="absolute left-[-8%] right-[-8%] top-1/2 block"
                  style={{
                    height: '2px',
                    background: 'black',
                    opacity: 0.95,
                    transform: 'translateY(-90%)',
                  }}
                />
              </span>
            </span>
          );
        }

        // default letter with slight optical kerning on edges
        return (
          <span
            key={`${ch}-${i}`}
            className="inline-block"
            style={{
              paddingLeft: i === 0 ? '0.012em' : '0',
              paddingRight: i === chars.length - 1 ? '0.006em' : '0',
            }}
          >
            {ch}
          </span>
        );
      })}
    </div>
  );
}
