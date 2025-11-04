// components/chat/ChatUI.tsx - Enhanced with SSE and dynamic updates
'use client';
import { useEffect, useRef, useState, useCallback } from 'react';
import { useModel } from '@/components/context/ModelContext';
import { useJobs } from '@/components/context/JobsContext';
import { useChat, Role } from '@/components/context/ChatContext';
import { getSessionId } from '@/components/utils/session';
import { cn } from '@/components/utils/cn';
import AnimatedOrb from './AnimatedOrb';

const MODEL_GRADIENT: Record<string, string> = {
    'gemini-2.5-pro': 'bg-gradient-to-br from-sky-400 to-blue-600',
    'gemini-2.5-flash': 'bg-gradient-to-br from-fuchsia-400 to-purple-600',
    'gemini-2.5-flash-lite': 'bg-gradient-to-br from-emerald-400 to-teal-600',
    'local': 'bg-gradient-to-br from-amber-400 to-orange-600',
};

const BACKEND_URL = process.env.NEXT_PUBLIC_BACKEND_URL;

interface StreamEvent {
    type: 'metadata' | 'message' | 'status' | 'error' | 'summary' | 'done';
    content?: string;
    step?: string;
    node?: string;
    error?: string | null;
    session_id?: string;
    intent?: string;
    state?: any;
}

interface MessageWithDetails {
    role: Role;
    content: string;
    details?: StreamEvent[];
    isThinking?: boolean;
    currentStep?: string;
}

export default function ChatUI() {
    const { model } = useModel();
    const { addJob, updateJob } = useJobs();
    const { conversations, active, activeId, createConversation, appendMessage } = useChat();
    const sessionId = getSessionId();
    const [pendingJobQuery, setPendingJobQuery] = useState<string | null>(null);

    // Initialize conversation
    const bootRef = useRef(false);
    useEffect(() => {
        if (bootRef.current) return;
        if (conversations.length === 0) {
            bootRef.current = true;
            createConversation({
                model,
                greeting: 'Hi! I\'m F.A.D.E, your AI drug discovery assistant. Ask me about proteins, diseases, or describe the drug you\'d like to design!'
            });
        }
    }, [conversations.length, createConversation, model]);

    const [input, setInput] = useState('');
    const [isThinking, setIsThinking] = useState(false);
    const [currentThinkingMessage, setCurrentThinkingMessage] = useState<MessageWithDetails | null>(null);
    const listRef = useRef<HTMLDivElement>(null);
    const eventSourceRef = useRef<EventSource | null>(null);

    // Auto-scroll
    useEffect(() => {
        listRef.current?.scrollTo({ top: listRef.current.scrollHeight, behavior: 'smooth' });
    }, [active?.messages?.length, isThinking, currentThinkingMessage?.currentStep]);

    // Cleanup SSE on unmount
    useEffect(() => {
        return () => {
            if (eventSourceRef.current) {
                eventSourceRef.current.close();
            }
        };
    }, []);

    function handleNewConversation() {
        setInput('');
        setIsThinking(false);
        setCurrentThinkingMessage(null);
        setPendingJobQuery(null);
        if (eventSourceRef.current) {
            eventSourceRef.current.close();
            eventSourceRef.current = null;
        }
        createConversation({
            model,
            greeting: 'Hi! I\'m F.A.D.E, your AI drug discovery assistant. What can I help you discover today?'
        });
    }

    async function handleSSERequest(query: string) {
        if (!BACKEND_URL) {
            appendMessage({
                role: 'assistant',
                content: '❌ Backend not configured. Please set NEXT_PUBLIC_BACKEND_URL in your environment.'
            });
            return;
        }

        // Close any existing SSE connection
        if (eventSourceRef.current) {
            eventSourceRef.current.close();
        }

        // Add user message
        appendMessage({ role: 'user', content: query });
        
        // Create thinking message
        const thinkingMsg: MessageWithDetails = {
            role: 'assistant',
            content: '',
            isThinking: true,
            currentStep: 'Initializing...',
            details: []
        };
        setCurrentThinkingMessage(thinkingMsg);
        setIsThinking(true);

        try {
            const requestBody = {
                query: query,
                user_id: sessionId,
                session_id: sessionId  // Add session_id for proper session tracking
            };

            // Use fetch with SSE parsing
            const response = await fetch(`${BACKEND_URL.replace(/\/$/, '')}/chat-stream`, {
                method: 'POST',
                headers: { 
                    'Content-Type': 'application/json',
                },
                body: JSON.stringify(requestBody),
            });

            if (!response.ok) {
                throw new Error(`API error: ${response.status}`);
            }

            const reader = response.body?.getReader();
            const decoder = new TextDecoder();
            let buffer = '';

            if (!reader) {
                throw new Error('No response body');
            }

            while (true) {
                const { done, value } = await reader.read();
                if (done) break;

                buffer += decoder.decode(value, { stream: true });
                const lines = buffer.split('\n');
                buffer = lines.pop() || '';

                for (const line of lines) {
                    if (line.startsWith('event:')) {
                        // Parse SSE event
                        const eventType = line.slice(6).trim();
                        const nextLine = lines[lines.indexOf(line) + 1];
                        
                        if (nextLine?.startsWith('data:')) {
                            const dataStr = nextLine.slice(5).trim();
                            
                            if (dataStr && dataStr !== '{}') {
                                try {
                                    const data = JSON.parse(dataStr);
                                    handleStreamEvent(eventType, data, thinkingMsg);
                                } catch (e) {
                                    console.error('Failed to parse SSE data:', e, dataStr);
                                }
                            }
                        }
                    }
                }
            }
        } catch (error) {
            console.error('SSE error:', error);
            setIsThinking(false);
            setCurrentThinkingMessage(null);
            appendMessage({
                role: 'assistant',
                content: `❌ Connection error: ${error instanceof Error ? error.message : 'Unknown error'}`
            });
        }
    }

    function handleStreamEvent(eventType: string, data: any, thinkingMsg: MessageWithDetails) {
        // Store all events in details for "Show Your Work"
        const event: StreamEvent = { type: eventType as any, ...data };
        thinkingMsg.details?.push(event);

        switch (eventType) {
            case 'metadata':
                // Session info, store but don't display
                console.log('Session:', data.session_id, 'Intent:', data.intent);
                break;

            case 'message':
                // Update current step display
                if (data.step) {
                    thinkingMsg.currentStep = data.content || `Processing ${data.step}...`;
                    setCurrentThinkingMessage({ ...thinkingMsg });
                }
                break;

            case 'status':
                // Update step based on status
                if (data.step) {
                    const stepMessages: Record<string, string> = {
                        'research_complete': '✅ Research complete',
                        'validation_complete': '✅ Target validated',
                        'structure_found': '✅ Structure resolved',
                        'pockets_detected': '✅ Pockets identified',
                        'molecules_generated': '✅ Molecules generated',
                        'screening_complete': '✅ Screening finished',
                        'pocket_detection_failed': '⚠️ Pocket detection failed',
                    };
                    thinkingMsg.currentStep = stepMessages[data.step] || data.step;
                    setCurrentThinkingMessage({ ...thinkingMsg });
                }
                break;

            case 'error':
                // Show error in thinking bubble
                thinkingMsg.currentStep = `❌ Error: ${data.error}`;
                setCurrentThinkingMessage({ ...thinkingMsg });
                break;

            case 'summary':
                // Final summary - replace thinking bubble with final message
                setIsThinking(false);
                setCurrentThinkingMessage(null);
                
                // Add final message with all details stored
                appendMessage({
                    role: 'assistant',
                    content: data.content,
                    details: thinkingMsg.details
                } as any);

                // If it's a job, add to jobs context
                if (data.job_id) {
                    addJob({
                        id: data.job_id,
                        query: input,
                        model: model,
                        user_id: sessionId,
                        current_time: new Date().toISOString(),
                        status: 'running'
                    });
                }
                break;

            case 'done':
                // Stream complete
                if (isThinking) {
                    // If still thinking and we get done, show what we have
                    setIsThinking(false);
                    setCurrentThinkingMessage(null);
                    
                    if (thinkingMsg.details && thinkingMsg.details.length > 0) {
                        appendMessage({
                            role: 'assistant',
                            content: thinkingMsg.currentStep || 'Process completed.',
                            details: thinkingMsg.details
                        } as any);
                    }
                }
                break;
        }
    }

    async function onSend() {
        if (!input.trim() || isThinking) return;

        const query = input.trim();
        setInput('');

        // Use SSE for all queries
        await handleSSERequest(query);
    }

    return (
        <div className="relative space-y-4">
            <div className="h-[60vh] w-full space-y-2 overflow-y-auto rounded-xl" ref={listRef}>
                {!active || active.messages.length === 0 ? (
                    <>
                        <AnimatedOrb gradient={MODEL_GRADIENT[model] || MODEL_GRADIENT.local} />
                        <div className="space-y-3 px-1">
                            <div className="text-sm text-white/70">
                                Try one of these drug discovery queries:
                            </div>
                            <QuickActionButton
                                onClick={async () => {
                                    const query = "Find EGFR inhibitors for lung cancer";
                                    setInput(query);
                                    await handleSSERequest(query);
                                }}
                                disabled={isThinking}
                            >
                                🎯 Find EGFR inhibitors for lung cancer
                            </QuickActionButton>
                            <QuickActionButton
                                onClick={async () => {
                                    const query = "What is KRAS and why is it important?";
                                    setInput(query);
                                    await handleSSERequest(query);
                                }}
                                disabled={isThinking}
                            >
                                🧬 What is KRAS and why is it important?
                            </QuickActionButton>
                            <QuickActionButton
                                onClick={async () => {
                                    const query = "Find PDB structure 6DUK";
                                    setInput(query);
                                    await handleSSERequest(query);
                                }}
                                disabled={isThinking}
                            >
                                🔬 Find PDB structure 6DUK
                            </QuickActionButton>
                        </div>
                    </>
                ) : (
                    <>
                        {active.messages.map((msg, i) => (
                            <EnhancedBubble key={i} message={msg as MessageWithDetails} />
                        ))}
                        {currentThinkingMessage && (
                            <ThinkingBubble 
                                currentStep={currentThinkingMessage.currentStep || 'Processing...'} 
                                details={currentThinkingMessage.details}
                            />
                        )}
                    </>
                )}
            </div>

            <ChatBar
                value={input}
                onInput={setInput}
                onSend={onSend}
                disabled={isThinking}
                placeholder={isThinking ? "F.A.D.E is thinking..." : "Ask about proteins, drugs, or diseases..."}
            />

            <StatusBar
                count={conversations.length}
                modelLabel={model}
                onNew={handleNewConversation}
            />
        </div>
    );
}

// New component for the dynamic thinking bubble
function ThinkingBubble({ currentStep, details }: { currentStep: string; details?: StreamEvent[] }) {
    return (
        <div className="flex justify-start">
            <div className="max-w-[85%] rounded-2xl bg-white/7 px-4 py-3 text-sm leading-relaxed text-white border border-white/10">
                <div className="flex items-center gap-3">
                    <div className="flex gap-1">
                        <span className="inline-block h-2 w-2 animate-pulse rounded-full bg-blue-400" />
                        <span className="inline-block h-2 w-2 animate-pulse rounded-full bg-blue-400 animation-delay-150" />
                        <span className="inline-block h-2 w-2 animate-pulse rounded-full bg-blue-400 animation-delay-300" />
                    </div>
                    <span className="text-white/90">{currentStep}</span>
                </div>
            </div>
        </div>
    );
}

// Enhanced bubble with "Show Your Work" toggle
function EnhancedBubble({ message }: { message: MessageWithDetails }) {
    const [showDetails, setShowDetails] = useState(false);
    const hasDetails = message.details && message.details.length > 0;

    if (message.role === 'user') {
        return <Bubble role={message.role} text={message.content} />;
    }

    return (
        <div className="flex justify-start">
            <div className="max-w-[85%] space-y-2">
                <div className="rounded-2xl bg-white/7 px-4 py-3 text-sm leading-relaxed text-white border border-white/10">
                    <MarkdownText text={message.content} />
                    
                    {hasDetails && (
                        <button
                            onClick={() => setShowDetails(!showDetails)}
                            className="mt-3 flex items-center gap-2 text-xs text-white/60 hover:text-white/80 transition-colors"
                        >
                            <svg 
                                className={cn("h-3 w-3 transition-transform", showDetails && "rotate-90")}
                                viewBox="0 0 24 24" 
                                fill="none" 
                                stroke="currentColor" 
                                strokeWidth="2"
                            >
                                <path d="M9 18l6-6-6-6" />
                            </svg>
                            {showDetails ? 'Hide' : 'Show'} Details
                        </button>
                    )}
                </div>

                {showDetails && hasDetails && (
                    <div className="rounded-lg bg-black/30 border border-white/5 p-3 text-xs font-mono text-white/70 max-h-64 overflow-y-auto">
                        {message.details!.map((event, i) => (
                            <div key={i} className="mb-1">
                                {event.type === 'message' && (
                                    <span className="text-blue-400">[message]</span>
                                )}
                                {event.type === 'status' && (
                                    <span className="text-green-400">[status]</span>
                                )}
                                {event.type === 'error' && (
                                    <span className="text-red-400">[error]</span>
                                )}
                                {' '}
                                {event.content || event.error || `${event.node || ''} ${event.step || ''}`.trim()}
                            </div>
                        ))}
                    </div>
                )}
            </div>
        </div>
    );
}

// Rest of the helper components remain the same
function ChatBar({ 
    value, 
    onInput, 
    onSend, 
    disabled,
    placeholder 
}: { 
    value: string;
    onInput: (v: string) => void;
    onSend: () => void;
    disabled?: boolean;
    placeholder?: string;
}) {
    return (
        <div className="relative flex items-center gap-2 rounded-full border border-white/10 bg-white/5 px-5 py-3 transition-colors focus-within:border-white/20 focus-within:bg-white/7">
            <input
                type="text"
                value={value}
                onChange={(e) => onInput(e.target.value)}
                onKeyDown={(e) => {
                    if (e.key === 'Enter' && !e.shiftKey && !disabled) {
                        e.preventDefault();
                        onSend();
                    }
                }}
                placeholder={placeholder}
                disabled={disabled}
                className="flex-1 bg-transparent text-sm text-white placeholder-white/40 outline-none disabled:opacity-50"
                aria-label="Message input"
            />
            <button
                type="button"
                onClick={onSend}
                disabled={disabled || !value.trim()}
                className="rounded-full bg-blue-500 p-2 text-white hover:bg-blue-600 disabled:opacity-50 disabled:cursor-not-allowed transition-colors"
                aria-label="Send message"
            >
                <svg viewBox="0 0 24 24" className="h-4 w-4" fill="none" stroke="currentColor" strokeWidth="1.8">
                    <path d="M22 2L11 13M22 2l-7 20-4-9-9-4 20-7z" />
                </svg>
            </button>
        </div>
    );
}

function StatusBar({ count, modelLabel, onNew }: { count: number; modelLabel: string; onNew: () => void }) {
    return (
        <div className="flex items-center justify-between text-xs text-white/60">
            <div>
                Conversations: <span className="text-white">{count}</span> · Engine: <span className="text-white">{modelLabel}</span>
            </div>
            <button
                type="button"
                onClick={onNew}
                className="inline-flex items-center gap-2 rounded-full border border-white/15 bg-white/10 px-3 py-1.5 text-xs font-medium text-white hover:bg-white/15 transition-colors"
                aria-label="Start new conversation"
                title="Start new conversation"
            >
                <svg viewBox="0 0 24 24" className="h-3.5 w-3.5" fill="none" stroke="currentColor" strokeWidth="1.8">
                    <path d="M12 5v14M5 12h14" />
                </svg>
                New conversation
            </button>
        </div>
    );
}

function Bubble({ role, text }: { role: Role; text: string }) {
    const isUser = role === 'user';

    return (
        <div className={cn('flex w-full', isUser ? 'justify-end' : 'justify-start')}>
            <div className={cn(
                'max-w-[85%] rounded-2xl px-4 py-3 text-sm leading-relaxed',
                isUser
                    ? 'bg-blue-500 text-white'
                    : 'bg-white/7 text-white border border-white/10'
            )}>
                <span>{text}</span>
            </div>
        </div>
    );
}

function MarkdownText({ text }: { text: string }) {
    // Simple markdown parsing for **bold** and line breaks
    const parts = text.split(/(\*\*[^*]+\*\*)/g);

    return (
        <div className="space-y-2">
            {parts.map((part, index) => {
                if (part.startsWith('**') && part.endsWith('**')) {
                    return (
                        <span key={index} className="font-semibold">
                            {part.slice(2, -2)}
                        </span>
                    );
                } else {
                    return part.split('\n').map((line, lineIndex) => (
                        <div key={`${index}-${lineIndex}`}>
                            {line}
                        </div>
                    ));
                }
            })}
        </div>
    );
}

function QuickActionButton({
    children,
    onClick,
    disabled
}: {
    children: React.ReactNode;
    onClick: () => void;
    disabled?: boolean;
}) {
    return (
        <button
            onClick={onClick}
            disabled={disabled}
            className="w-full text-left rounded-xl border border-white/10 bg-white/5 px-4 py-3 text-sm text-white/90 hover:bg-white/10 hover:border-white/20 disabled:opacity-50 disabled:cursor-not-allowed transition-colors"
        >
            {children}
        </button>
    );
}
