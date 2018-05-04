/* 
 * File:   Thread.h
 * Author: Null
 * Blog: http://hi.baidu.com/hetaoos
 * Created on 2008��7��30��, ����10:13
 */

/*
 * �ڱ����ʱ��ǵü��ϲ�����-lpthread
 * 
 */

#ifndef _THREAD_H
#define _THREAD_H

#include <pthread.h>
#include <unistd.h>
/*
 * �߳�����ʵ����
 */
class Runnable
{
	public:
		//����ʵ��
		virtual void run() = 0;
};

/*
 * �߳���
 */
class Thread : public Runnable
{
	private:
		//�̳߳�ʼ�����
		static int threadInitNumber;
		//��ǰ�̳߳�ʼ�����
		int curThreadInitNumber;
		//�߳���
		Runnable *target;
		//��ǰ�̵߳��߳�ID
		pthread_t tid;
		//�̵߳�״̬
		int threadStatus;
		//�߳�����
		pthread_attr_t attr;
		//�߳����ȼ�
		sched_param param;
		//��ȡִ�з�����ָ��
		static void* run0(void* pVoid);
		//�ڲ�ִ�з���
		void* run1();
		//��ȡһ���߳����
		static int getNextThreadNum();

	public:
		//�̵߳�״̬���½�
		static const int THREAD_STATUS_NEW = 0;
		//�̵߳�״̬����������
		static const int THREAD_STATUS_RUNNING = 1;
		//�̵߳�״̬�����н���
		static const int THREAD_STATUS_EXIT = -1;
		//���캯��
		Thread();
		//���캯��
		Thread(Runnable *iTarget);
		//����
		~Thread();
		//�̵߳�����ʵ��
		void run();
		//��ʼִ���߳�
		bool start();
		//��ȡ�߳�״̬
		int getState();
		//�ȴ��߳�ֱ���˳�
		void join();
		//�ȴ��߳��˳����߳�ʱ
		void join(unsigned long millisTime);
		//�Ƚ������߳�ʱ����ͬ,ͨ�� curThreadInitNumber �ж�
		bool operator ==(const Thread *otherThread);
		//��ȡThis�߳�ID
		pthread_t getThreadID();
		//��ȡ��ǰ�߳�ID
		static pthread_t getCurrentThreadID();
		//��ǰ�߳��Ƿ��ĳ���߳����,ͨ�� tid �ж�
		static bool isEquals(Thread *iTarget);
		//�����̵߳�����:��/�ǰ�
		void setThreadScope(bool isSystem);
		//��ȡ�̵߳�����:��/�ǰ�
		bool getThreadScope();
		//�����̵߳����ȼ�,1-99,����99Ϊʵʱ;�����Ϊ��ͨ
		void setThreadPriority(int priority);
		//��ȡ�̵߳����ȼ�
		int getThreadPriority();
};


int Thread::threadInitNumber = 1;

int Thread::getNextThreadNum()
{
	return threadInitNumber++;
}

void* Thread::run0(void* pVoid)
{
	Thread* p = (Thread*) pVoid;
	p->run1();
	return p;
}

void* Thread::run1()
{

	threadStatus = THREAD_STATUS_RUNNING;
	tid = pthread_self();
	run();
	threadStatus = THREAD_STATUS_EXIT;
	tid = 0;
	pthread_exit(NULL);
}

void Thread::run()
{
	if (target != NULL)
	{
		(*target).run();
	}
}

Thread::Thread()
{
	tid = 0;
	threadStatus = THREAD_STATUS_NEW;
	curThreadInitNumber = getNextThreadNum();
	pthread_attr_init(&attr);
}

Thread::Thread(Runnable *iTarget)
{
	target = iTarget;
	tid = 0;
	threadStatus = THREAD_STATUS_NEW;
	curThreadInitNumber = getNextThreadNum();
	pthread_attr_init(&attr);
}

Thread::~Thread()
{
	pthread_attr_destroy(&attr);
}

bool Thread::start()
{
	return pthread_create(&tid, &attr, run0, this) == 0;
}

pthread_t Thread::getCurrentThreadID()
{
	return pthread_self();
}

pthread_t Thread::getThreadID()
{
	return tid;
}

int Thread::getState()
{
	return threadStatus;
}

void Thread::join()
{
	if (tid > 0)
	{
		pthread_join(tid, NULL);
	}
}

void Thread::join(unsigned long millisTime)
{

	if (tid == 0)
	{
		return;
	}
	if (millisTime == 0)
	{
		join();
	}else
	{
		unsigned long k = 0;
		while (threadStatus != THREAD_STATUS_EXIT && k <= millisTime)
		{
			usleep(100);
			k++;
		}
	}
}

bool Thread::operator ==(const Thread *otherThread)
{
	if (otherThread == NULL)
	{
		return false;
	}
	if (curThreadInitNumber == (*otherThread).curThreadInitNumber)
	{
		return true;
	}
	return false;
}

bool Thread::isEquals(Thread *iTarget)
{
	if (iTarget == NULL)
	{
		return false;
	}
	return pthread_self() == iTarget->tid;
}

void Thread::setThreadScope(bool isSystem)
{
	if (isSystem)
	{
		pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);
	}else
	{
		pthread_attr_setscope(&attr, PTHREAD_SCOPE_PROCESS);
	}
}

bool Thread::getThreadScope()
{
	int scopeType = 0;
	pthread_attr_getscope(&attr, &scopeType);
	return scopeType == PTHREAD_SCOPE_SYSTEM;
}

void Thread::setThreadPriority(int priority)
{
	pthread_attr_getschedparam(&attr, &param);
	param.__sched_priority = priority;
	pthread_attr_setschedparam(&attr, &param);
}

int Thread::getThreadPriority()
{
	pthread_attr_getschedparam(&attr, &param);
	return param.__sched_priority;
}
#endif /* _THREAD_H */
