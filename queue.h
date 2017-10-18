#ifndef _QUEUE_H_
#define _QUEUE_H_

#include "debug.h"
#include <cstddef>

template <typename T>
class queue {
	T *data;
	ssize_t head;
	ssize_t tail;
	ssize_t size;
	void grow();
	void check() const;
public:
	queue();
	queue(const queue &);
	~queue();
	ssize_t length() const;
	bool empty() const;
	void enqueue(const T &);
	void dequeue(T &);
	T const dequeue();
	queue const &operator=(const queue &);
};

template <typename T>
queue<T>::queue() : data(0), head(-1), tail(-1), size(0)
{
}

template <typename T>
queue<T>::queue(queue const &q) : data(0), head(-1), tail(-1), size(0)
{
	if (q.empty())
		return;

	ssize_t n = 0;
	ssize_t h = q.head;
	data = new T[q.length() + 1];
	while (h != q.tail) {
		ASSERT(n < q.length());
		data[n] = q.data[h];
		h = (h + 1) % q.size;
		n++;
	}
	ASSERT(n == q.length());
	head = 0;
	tail = q.length();
	size = q.length() + 1;
	check();
}

template <typename T>
queue<T>::~queue()
{
	check();
	delete[] data;
}

template <typename T>
ssize_t queue<T>::length() const
{
	check();
	if (size == 0)
		return 0;
	return (tail >= head) ? tail - head
	                      : size + tail - head;
}

template <typename T>
bool queue<T>::empty() const
{
	check();
	return (head == tail) ? true : false;
}

template <typename T>
void queue<T>::enqueue(T const &t)
{
	check();
	if (length() + 1 >= size)
		grow();
	data[tail] = t;
	tail = (tail + 1) % size;
	check();
}

template <typename T>
void queue<T>::dequeue(T &t)
{
	check();
	if (empty())
		std::abort();
	t = data[head];
	head = (head + 1) % size;
	check();
}

template <typename T>
T const queue<T>::dequeue()
{
	ssize_t what = head;
	check();
	if (empty())
		std::abort();
	head = (head + 1) % size;
	check();
	return data[what];
}

template <typename T>
queue<T> const &queue<T>::operator=(const queue &q)
{
	check();
	q.check();

	if (this != &q) {
		ssize_t n = 0;
		ssize_t h = q.head;

		if (size <= q.length()) {
			delete[] data;
			data = new T[size = q.length() + 1];
		}

		while (h != q.tail) {
			data[n] = q.data[h];
			h = (h + 1) % q.size;
			n++;
		}

		head = 0;
		tail = q.length();
	}

	check();
	ASSERT(length() == q.length());

	return *this;
}

template <typename T>
void queue<T>::grow() 
{
	if (!data) {
		data = new T[2];
		head = 0;
		tail = 0;
		size = 2;
	} else {
		ssize_t n = 0;
		ssize_t h = head;
		ssize_t s = 2 * size;
		T *d = new T[s];
		while (h != tail) {
			ASSERT(n < s);
			d[n] = data[h];
			h = (h + 1) % size;
			n++;
		}
		ASSERT(n == length());
		delete[] data, data = d;
		head = 0;
		tail = n;
		size = s;
	}
	check();
}

template <typename T>
void queue<T>::check() const
{
	ASSERT((!data
	        && !size)
	       || (data 
	           && size > 0 
	           && head >= 0 
	           && head < size
	           && tail >= 0
	           && tail < size));
}

#endif
