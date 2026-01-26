#ifndef THREAD_POOL_HPP
#define THREAD_POOL_HPP

#include "std_include.hpp"

namespace thread_pool{

using std::size_t;
using std::views::iota;
using std::vector;

template<typename T> concept NonVoid = !std::same_as<T, void>;

template<NonVoid score_t> struct Instruction{
	vector<std::variant<std::function<void(size_t)>, std::function<score_t(size_t)> > > mapFuncs;
	vector<std::function<score_t(score_t, score_t)> > reduceFuncs;
	vector<score_t> zeros;
};

template<NonVoid score_t, std::ranges::range Range> struct Task{
	Instruction<score_t> const *const instr;
	Range elementRange;
	
	Task(Instruction<score_t> const *const instr, Range elementRange) noexcept:
		instr(instr), elementRange(elementRange){}
	
	vector<score_t> operator()() const noexcept{
		vector<score_t> res;
		size_t iReduce = 0;
		for (size_t iMap : iota((size_t) 0, instr->mapFuncs.size())){
			auto map = instr->mapFuncs[iMap];
			if (std::holds_alternative<std::function<void(size_t)> >(map)){
				for (size_t iElement: elementRange){
					std::get<std::function<void(size_t)> >(map)(iElement);
				}
			}
			else{
				auto reduce = instr->reduceFuncs[iReduce];
				score_t temp = instr->zeros[iReduce];
				for (size_t iElement: elementRange){
					temp = reduce(temp, std::get<std::function<score_t(size_t)> >(map)(iElement));
				}
				res.push_back(temp);
				iReduce++;
			}
		}
		return res;
	}
};


template<template<typename> class Scheduler, typename score_t> concept SCHEDULER = requires(Scheduler<score_t> scheduler, Instruction<score_t> const *const instr, size_t index){
    {Scheduler<score_t>{instr, index, index, index}} noexcept;
    {scheduler(index)} noexcept;
	{scheduler.getResults()} noexcept -> std::same_as<vector<vector<score_t> > >;
	requires NonVoid<score_t>;
};

template<typename score_t> class SimpleScheduler{
	vector<Task<score_t, std::ranges::iota_view<size_t, size_t> > > tasks;
	vector<vector<score_t> > results;

public:
	SimpleScheduler(Instruction<score_t> const *const instr, size_t nThreads, size_t iElementBegin, size_t iElementEnd) noexcept: results(nThreads){
		for (size_t iThread : iota((size_t) 0, nThreads)){
			tasks.emplace_back(instr, iota(iElementBegin + iThread * (iElementEnd - iElementBegin) / nThreads, iElementBegin + (iThread + 1) * (iElementEnd - iElementBegin) / nThreads));
		}
	}
	
	void operator()(size_t iThread) noexcept{
		results[iThread] = tasks[iThread]();
	}
	
	vector<vector<score_t> > getResults() noexcept{
		return std::move(results);
	}
};

template<typename TP> concept THREAD_POOL = requires(TP tp, Instruction<typename TP::score_t> const &instr, size_t index){
	{ TP{ index, index, index } } noexcept;
	{tp(instr)} noexcept -> std::same_as<vector<typename TP::score_t> >;
};

template<typename T, template<typename> class Scheduler> requires SCHEDULER<Scheduler, T> class ThreadPool{
public:
	using score_t = T;

private:
	struct Block{
		std::unique_ptr<Scheduler<score_t> > scheduler;
		std::atomic<size_t> readyCnt;
		std::atomic<bool> dataReady;
		std::shared_ptr<Block> next;
	};
	
	static void work(std::shared_ptr<Block> block, size_t iThread) noexcept{
		while(block.get() != nullptr){
			while(!block->dataReady);
			if (block->scheduler.get() != nullptr) (*block->scheduler)(iThread);
			block->readyCnt++;
			block = block->next;
		}
	}
	
	std::shared_ptr<Block> block;
	vector<std::thread> threads;
	size_t nThreads, iElementBegin, iElementEnd;

public:
	ThreadPool(size_t nThreads, size_t iElementBegin, size_t iElementEnd) noexcept: nThreads(nThreads), iElementBegin(iElementBegin), iElementEnd(iElementEnd), block(new Block){
		for (size_t iThread : iota((size_t) 1, nThreads)){
			threads.emplace_back(&ThreadPool<score_t, Scheduler>::work, block, iThread);
		}
	}
	
	~ThreadPool() noexcept{
		block->dataReady = true;
		for (std::thread &thread : threads){
			thread.join();
		}
	}
	
	vector<score_t> operator()(Instruction<score_t> const &instr) noexcept{
		if (nThreads == 1) return Task<score_t, std::ranges::iota_view<size_t, size_t> >(&instr, iota(iElementBegin, iElementEnd))();
		block->scheduler.reset(new Scheduler<score_t>(&instr, nThreads, iElementBegin, iElementEnd));
		block->next.reset(new Block);
		block->dataReady = true;
		if (block->scheduler.get() != nullptr) (*block->scheduler)(0);
		std::unique_lock<std::mutex> lk(block->headMtx);
		while (block->readyCnt < nThreads - 1);
		vector<vector<score_t> > results = block->scheduler->getResults();
		vector<score_t> res;
		for (size_t iReduce : iota((size_t) 0, instr.reduceFuncs.size())){
			auto reduce = instr.reduceFuncs[iReduce];
			score_t temp = instr.zeros[iReduce];
			for (size_t iThread: iota((size_t) 0, nThreads)){
				temp = reduce(temp, results[iThread][iReduce]);
			}
			res.push_back(temp);
		}
		block = block->next;
		return res;
	}
};

};
#endif