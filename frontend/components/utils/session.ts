// components/utils/session.ts
export function getSessionId(key = 'fade_session_id') {
  if (typeof window === 'undefined') return 'server';
  let sid = localStorage.getItem(key);
  if (!sid) {
    sid = crypto.randomUUID();
    localStorage.setItem(key, sid);
  }
  return sid;
}
