// components/chat/ChatUI.tsx - Updated to connect with F.A.D.E Gateway
'use client';
import { useEffect, useRef, useState } from 'react';
import { useModel } from '@/components/context/ModelContext';
import { useJobs } from '@/components/context/JobsContext';
import { useChat, Role } from '@/components/context/ChatContext';
import { getSessionId } from '@/components/utils/session';
import { cn } from '@/components/utils/cn';
import AnimatedOrb from './AnimatedOrb';

const MODEL_GRADIENT: Record<string, string> = {
    'gpt-4o': 'bg-gradient-to-br from-sky-400 to-blue-600',
    'gpt-4.1': 'bg-gradient-to-br from-fuchsia-400 to-purple-600',
    'gpt-3.5': 'bg-gradient-to-br from-emerald-400 to-teal-600',
    'local': 'bg-gradient-to-br from-amber-400 to-orange-600',
};

const BACKEND_URL = process.env.NEXT_PUBLIC_BACKEND_URL;

interface GatewayResponse {
    type: 'chat' | 'job' | 'error';
    message?: string;
    job_id?: string;
    status?: string;
    reasoning?: string;
    confidence?: number;
    error?: string;
}

export default function ChatUI() {
    const { model } = useModel();
    const { addJob, updateJob } = useJobs();
    const { conversations, active, activeId, createConversation, appendMessage } = useChat();
    const sessionId = getSessionId();

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
    const listRef = useRef<HTMLDivElement>(null);

    // Auto-scroll
    useEffect(() => {
        listRef.current?.scrollTo({ top: listRef.current.scrollHeight, behavior: 'smooth' });
    }, [active?.messages.length, isThinking]);

    function handleNewConversation() {
        setInput('');
        setIsThinking(false);
        createConversation({ 
            model, 
            greeting: 'Hi! I\'m F.A.D.E, your AI drug discovery assistant. What can I help you discover today?' 
        });
    }

    async function handleGatewayRequest(query: string) {
        if (!BACKEND_URL) {
            appendMessage({ 
                role: 'assistant', 
                content: '❌ Gateway not configured. Please set NEXT_PUBLIC_BACKEND_URL in your environment.' 
            });
            return;
        }

        const requestPayload = {
            id: crypto.randomUUID(),
            query: query,
            model: model,
            current_time: new Date().toISOString(),
            user_id: sessionId
        };

        try {
            console.log('📤 Sending to gateway:', requestPayload);

            const response = await fetch(`${BACKEND_URL.replace(/\/$/, '')}/jobs`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify(requestPayload),
            });

            if (!response.ok) {
                throw new Error(`Gateway error: ${response.status} ${response.statusText}`);
            }

            const result: GatewayResponse = await response.json();
            console.log('📥 Gateway response:', result);

            if (result.type === 'chat') {
                // Display chat response directly
                let message = result.message || 'I received your question but couldn\'t generate a response.';
                
                // Add confidence indicator if available
                if (result.confidence && result.confidence < 0.8) {
                    message += `\n\n*Note: I'm ${Math.round(result.confidence * 100)}% confident about this response.*`;
                }

                appendMessage({ 
                    role: 'assistant', 
                    content: message 
                });

            } else if (result.type === 'job') {
                // Handle job submission
                const jobId = result.job_id!;
                const jobData = {
                    ...requestPayload,
                    id: jobId,
                    status: result.status as any || 'queued'
                };

                // Add to job queue
                addJob(jobData);

                // Show job confirmation in chat
                let jobMessage = `🚀 **Drug Discovery Pipeline Started!**\n\n`;
                jobMessage += `**Job ID:** ${jobId.slice(0, 8)}...\n`;
                jobMessage += `**Status:** ${result.status}\n`;
                jobMessage += `**Query:** "${query}"\n\n`;
                jobMessage += `Your request has been submitted to the F.A.D.E pipeline. `;
                jobMessage += `You can monitor progress in the "In Progress" tab.\n\n`;
                
                if (result.confidence && result.confidence < 0.9) {
                    jobMessage += `*Classification confidence: ${Math.round(result.confidence * 100)}%*`;
                }

                appendMessage({ 
                    role: 'assistant', 
                    content: jobMessage 
                });

                // Start monitoring job status
                monitorJobStatus(jobId);

            } else if (result.type === 'error') {
                // Handle errors
                appendMessage({ 
                    role: 'assistant', 
                    content: `❌ Sorry, I encountered an error: ${result.message || result.error || 'Unknown error'}` 
                });
            }

        } catch (error) {
            console.error('❌ Gateway request failed:', error);
            appendMessage({ 
                role: 'assistant', 
                content: `❌ Sorry, I couldn't connect to the gateway: ${error instanceof Error ? error.message : 'Unknown error'}. Please try again.` 
            });
        }
    }

    async function monitorJobStatus(jobId: string) {
        // Monitor job status and update chat when completed
        const checkStatus = async () => {
            try {
                const response = await fetch(`${BACKEND_URL}/jobs/${jobId}/status`);
                if (response.ok) {
                    const status = await response.json();
                    updateJob(jobId, { status: status.status, error: status.error });

                    if (status.status === 'completed') {
                        appendMessage({
                            role: 'assistant',
                            content: `✅ **Pipeline Completed!**\n\nJob ${jobId.slice(0, 8)}... has finished successfully. Check the "In Progress" tab to view detailed results.`
                        });
                    } else if (status.status === 'failed') {
                        appendMessage({
                            role: 'assistant', 
                            content: `❌ **Pipeline Failed**\n\nJob ${jobId.slice(0, 8)}... encountered an error: ${status.error || 'Unknown error'}`
                        });
                    } else if (status.status === 'running') {
                        // Continue monitoring
                        setTimeout(checkStatus, 10000); // Check every 10 seconds
                    }
                }
            } catch (error) {
                console.error('Status check failed:', error);
            }
        };

        // Start monitoring after a brief delay
        setTimeout(checkStatus, 5000);
    }

    async function onSend() {
        const text = input.trim();
        if (!text || !active || isThinking) return;

        // Add user message
        appendMessage({ role: 'user', content: text });
        setInput('');
        setIsThinking(true);

        // Send to gateway
        await handleGatewayRequest(text);
        setIsThinking(false);
    }

    async function onQuickAction(text: string) {
        if (isThinking || !active) return;

        // Add user message
        appendMessage({ role: 'user', content: text });
        setIsThinking(true);

        // Send to gateway
        await handleGatewayRequest(text);
        setIsThinking(false);
    }

    const messages = active?.messages ?? [];

    return (
        <div className="flex flex-col items-center">
            {/* Hero Section */}
            <div className="flex flex-col items-center text-center">
                <div className="relative">
                    <AnimatedOrb gradient={MODEL_GRADIENT[model] ?? MODEL_GRADIENT['gpt-4o']} />
                    <div className="absolute -top-1 -right-1 w-6 h-6 bg-green-500 rounded-full flex items-center justify-center">
                        <span className="text-xs">🧬</span>
                    </div>
                </div>
                <h1 className="mt-6 text-3xl md:text-4xl font-semibold">F.A.D.E Drug Discovery</h1>
                <p className="mt-2 text-white/70 max-w-2xl text-sm">
                    Ask me about proteins, diseases, or describe your drug discovery goals. I'll provide information or run the full pipeline for you.
                </p>
            </div>

            {/* Conversation Bar */}
            <div className="mt-6 w-full max-w-3xl">
                <NewConversationBar
                    onNew={handleNewConversation}
                    count={conversations.length}
                    modelLabel="F.A.D.E"
                />
            </div>

            {/* Chat Panel */}
            <div className="mt-8 w-full max-w-3xl">
                <div 
                    ref={listRef} 
                    className="h-[44vh] md:h-[40vh] overflow-y-auto rounded-2xl border border-white/10 bg-white/5 p-4"
                >
                    <div className="space-y-3">
                        {messages.map((m) => (
                            <Bubble key={m.id} role={m.role} text={m.content} />
                        ))}
                        {isThinking && <AssistantTyping />}
                    </div>
                </div>

                {/* Input Area */}
                <div className="mt-4 rounded-2xl border border-white/10 bg-white/5 p-3">
                    <div className="flex items-start gap-3">
                        <textarea
                            value={input}
                            onChange={(e) => setInput(e.target.value)}
                            onKeyDown={(e) => {
                                if (e.key === 'Enter' && !e.shiftKey) {
                                    e.preventDefault();
                                    onSend();
                                }
                            }}
                            rows={1}
                            placeholder="Ask about proteins, request drug discovery, or get pipeline status..."
                            className="min-h-[44px] max-h-40 flex-1 resize-y bg-transparent p-2 outline-none placeholder:text-white/50"
                        />
                        <button
                            onClick={onSend}
                            disabled={!input.trim() || isThinking || !activeId}
                            className="rounded-full bg-blue-500 px-4 py-2 text-sm font-medium text-white hover:bg-blue-600 disabled:opacity-50 disabled:cursor-not-allowed transition-colors"
                        >
                            {isThinking ? 'Thinking...' : 'Send'}
                        </button>
                    </div>
                </div>

                {/* Enhanced Quick Actions */}
                <div className="mt-4 space-y-2">
                    <QuickActionButton 
                        onClick={() => onQuickAction('What is F.A.D.E and how does it work?')}
                        disabled={isThinking}
                    >
                        💡 What is F.A.D.E and how does it work?
                    </QuickActionButton>
                    
                    <QuickActionButton 
                        onClick={() => onQuickAction('What is KRAS and why is it important in cancer?')}
                        disabled={isThinking}
                    >
                        🧬 What is KRAS and why is it important in cancer?
                    </QuickActionButton>
                    
                    <QuickActionButton 
                        onClick={() => onQuickAction('Find EGFR inhibitors for lung cancer with oral bioavailability')}
                        disabled={isThinking}
                    >
                        🎯 Find EGFR inhibitors for lung cancer (oral)
                    </QuickActionButton>
                    
                    <QuickActionButton 
                        onClick={() => onQuickAction('Design KRAS G12C inhibitors with brain penetration for pancreatic cancer')}
                        disabled={isThinking}
                    >
                        🧠 Design KRAS G12C inhibitors with brain penetration
                    </QuickActionButton>
                </div>
            </div>
        </div>
    );
}

// Supporting Components

function NewConversationBar({
    onNew,
    count,
    modelLabel,
}: {
    onNew: () => void;
    count: number;
    modelLabel: string;
}) {
    return (
        <div className="flex items-center justify-between rounded-xl border border-white/10 bg-white/5 px-3 py-2">
            <div className="text-xs text-white/70">
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
    const isMarkdown = text.includes('**') || text.includes('\n');
    
    return (
        <div className={cn('flex w-full', isUser ? 'justify-end' : 'justify-start')}>
            <div className={cn(
                'max-w-[85%] rounded-2xl px-4 py-3 text-sm leading-relaxed',
                isUser 
                    ? 'bg-blue-500 text-white' 
                    : 'bg-white/7 text-white border border-white/10'
            )}>
                {isMarkdown ? (
                    <MarkdownText text={text} />
                ) : (
                    <span>{text}</span>
                )}
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

function AssistantTyping() {
    return (
        <div className="flex justify-start">
            <div className="flex items-center gap-2 rounded-2xl bg-white/7 px-4 py-3 border border-white/10">
                <span className="text-sm text-white/80">F.A.D.E is thinking</span>
                <div className="flex gap-1">
                    <span className="inline-block h-1.5 w-1.5 animate-bounce rounded-full bg-blue-400" />
                    <span className="inline-block h-1.5 w-1.5 animate-bounce rounded-full bg-blue-400 animation-delay-150" />
                    <span className="inline-block h-1.5 w-1.5 animate-bounce rounded-full bg-blue-400 animation-delay-300" />
                </div>
            </div>
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