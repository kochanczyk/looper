// This file is a part of Looper, an implementation of trajectory looping,
// distributed under the GNU General Public License, version 3.0 (see file 
// License.txt contained in the source code package or the copy online at
// http://www.gnu.org/licenses/gpl-3.0.txt).

#ifndef CTHULHU_THREAD_POOL_HPP
#define CTHULHU_THREAD_POOL_HPP

#include <vector>
#include <queue>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <functional>
#include <stdexcept>


class ThreadPool {
public:
  ThreadPool(size_t);
  ~ThreadPool();

  ThreadPool(ThreadPool const&) = delete;
  ThreadPool& operator=(ThreadPool const&) = delete;

  template<class F, class... Args>
  auto enqueue(F&& f, Args&&... args) 
      -> std::future<typename std::result_of<F(Args...)>::type>;

protected:
  std::vector<std::thread> workers;
  std::queue<std::function<void()>> tasks;
  
  std::mutex queue_mutex;
  std::condition_variable condition;
  bool stop;
};
 

inline ThreadPool::ThreadPool(size_t threads)
: stop{false}
{
  for (size_t i = 0;i<threads;++i) {
    workers.emplace_back( 
      [this] {
        while (true) {
          std::function<void()> task;
          {
            std::unique_lock<std::mutex> lock(this->queue_mutex);
            this->condition.wait(lock, [this]{   return this->stop 
                                              || not this->tasks.empty(); });
            if(this->stop and this->tasks.empty()) { return; }
            task = std::move(this->tasks.front());
            this->tasks.pop();
          }
          task();
        } // while
      } // lambda body
    );// emplace_back
  } // for each thread
}


template<class F, class... Args>
auto ThreadPool::enqueue(F&& f, Args&&... args) 
    -> std::future<typename std::result_of<F(Args...)>::type>
{
  using return_type = typename std::result_of<F(Args...)>::type;
  auto task = std::make_shared<std::packaged_task<return_type()>>
              (std::bind(std::forward<F>(f), std::forward<Args>(args)...));

  std::future<return_type> res = task->get_future();
  {
    std::unique_lock<std::mutex> lock(queue_mutex);
    if (stop) { 
      throw std::runtime_error("Parallelization error :"
                               "Cannot enqueue to a stopped thread pool."); 
    }
    tasks.emplace( [task](){ (*task)(); } );
  }
  condition.notify_one();
  return res;
}


inline ThreadPool::~ThreadPool()
{
  {
    std::unique_lock<std::mutex> lock(queue_mutex);
    stop = true;
  }
  condition.notify_all();
  for (auto &w: workers) { w.join(); }
}

#endif

