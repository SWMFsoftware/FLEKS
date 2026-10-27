#ifndef _PARALLEL_H_
#define _PARALLEL_H_

#include <algorithm>
#include <cstdlib>
#include <exception>
#include <thread>
#include <vector>

// Number of worker threads for the text/compression work. It can be set with
// the environment variable CONVERTER_NUM_THREADS. The default is capped,
// because a large shared node may report hundreds of hardware threads, and
// every thread holds its own output buffer.
inline size_t converter_num_threads() {
  static const size_t n = [] {
    constexpr size_t maxDefaultThreads = 16;
    size_t nThread = std::min<size_t>(std::thread::hardware_concurrency(),
                                      maxDefaultThreads);
    if (const char* env = std::getenv("CONVERTER_NUM_THREADS")) {
      const long v = std::atol(env);
      if (v > 0)
        nThread = static_cast<size_t>(v);
    }
    return std::max<size_t>(nThread, 1);
  }();
  return n;
}

// Run f(0), ..., f(n-1) on up to n threads and wait for all of them.
// An exception thrown by any f(t) is rethrown on the calling thread after all
// threads have joined. If a thread cannot be started, its work runs on the
// calling thread instead.
template <typename F> void run_parallel(size_t n, F&& f) {
  if (n == 0)
    return;

  std::vector<std::exception_ptr> errors(n);
  auto task = [&](size_t t) {
    try {
      f(t);
    } catch (...) {
      errors[t] = std::current_exception();
    }
  };

  std::vector<std::thread> threads;
  size_t nStarted = 1; // Task 0 always runs on the calling thread.
  try {
    threads.reserve(n - 1);
    for (; nStarted < n; ++nStarted)
      threads.emplace_back([&task, t = nStarted] { task(t); });
  } catch (...) {
    // Out of threads or memory: fall through and run the rest here.
  }

  task(0);
  for (size_t t = nStarted; t < n; ++t)
    task(t);

  for (auto& th : threads)
    th.join();

  for (const auto& error : errors)
    if (error)
      std::rethrow_exception(error);
}

#endif
