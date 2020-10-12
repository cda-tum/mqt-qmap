/*
 * This file is part of the JKQ QMAP library which is released under the MIT license.
 * See file README.md or go to https://iic.jku.at/eda/research/ibm_qx_mapping/ for more information.
 */

#include <set>
#include <queue>
#include <cassert>
#include <iostream>
#include <vector>

#ifndef UNIQUE_PRIORITY_QUEUE_H
#define UNIQUE_PRIORITY_QUEUE_H

#define UNUSED(x)       {(void) x;}

constexpr int    MAX_QUEUE_SIZE               = 6000000;
constexpr int    MAX_NODES_MARGIN             = 500000;
constexpr int    MAX_QUEUE_COPY_LENGTH        = 1000000;
constexpr double QUEUE_COPY_LENGTH_PERCENTAGE = 1./6;

template<class T>
struct do_nothing {
    void operator()(const T&){
        // intentionally left blank
    }
};

template<class T, class Container = std::vector<T>, class Compare = std::less<typename Container::value_type>>
class own_priority_queue : public std::priority_queue<T, Container, Compare> {
	public:
		Container& get_container() {
            return this->c;
        }
};

/**
 * Priority queue with unique (according to FuncCompare) elements of type T where the sorting is based on CostCompare.
 * If NDEBUG is *not* defined, there are some assertions that help catching errors in the provided comparision functions.
 */
template<class T, class CostCompare = std::greater<T>, class FuncCompare = std::less<T>, class CleanObsoleteElement = do_nothing<T>>
class unique_priority_queue
{
public:
    typedef typename own_priority_queue<T, std::vector<T>, CostCompare>::size_type size_type;

    /**
     * Return true if the element was inserted into the queue.
     * This happens if equivalent element is present or if the new element has a lower cost associated to it.
     * False is returned if no insertion into the queue took place.
     */
    bool push(const T& v)
    {
        const auto& insertion_pair = membership_.insert(v);
        if(insertion_pair.second)
        {
            queue_.push(v);
        }
        else if (CostCompare()(*(insertion_pair.first), v))
        {
            CleanObsoleteElement()(*(insertion_pair.first));
            membership_.erase(insertion_pair.first);
       
            const auto inserted = membership_.insert(v);
            assert(inserted.second); UNUSED(inserted);
            
            queue_ = own_priority_queue<T, std::vector<T>, CostCompare>();
            for(const auto& element : membership_) {
                queue_.push(element);
            }
            assert(queue_.size() == membership_.size());
            
            return true;
        }
        else 
        {
            CleanObsoleteElement()(v);
        }
        assert(queue_.size() == membership_.size());
        return insertion_pair.second;
    }

    void pop()
    {
        assert(!queue_.empty() && queue_.size() == membership_.size());

        const auto& top_element   = queue_.top();
        const auto  number_erased = membership_.erase(top_element);

        assert(number_erased == 1); UNUSED(number_erased);

        queue_.pop();
        assert(queue_.size() == membership_.size());
    }

    const T& top() const
    {
        assert(!queue_.empty());
        return queue_.top();
    }

    bool empty() const
    {
        assert(queue_.size() == membership_.size());
        return queue_.empty();
    }

    size_type size() const
    {
        return queue_.size();
    }

    std::vector<T>& get_container() 
    {
        return queue_.get_container();
    }

    void delete_queue() 
    {
        std::vector<T> &v = get_container();
        for(typename std::vector<T>::const_iterator it = v.begin(); it != v.end();
            it++) {
            CleanObsoleteElement()(*it);
        }
        queue_      = own_priority_queue<T, std::vector<T>, CostCompare>();
        membership_ = std::set<T, FuncCompare>();
    }

    // clears the queue until a certain length is reached
    void update() 
    {	
        T temp_queue[MAX_QUEUE_SIZE];
        unsigned int length = std::min((int)(queue_.size() * QUEUE_COPY_LENGTH_PERCENTAGE), MAX_QUEUE_COPY_LENGTH);
        
        for(unsigned int i = 0; i < length; i++) {
            temp_queue[i] = queue_.top();
            queue_.pop();	

            last_node_copied = i;	
        }
        delete_queue();
        for(unsigned int i = 0; i < length; i++) {
            queue_.push(temp_queue[i]);		
        }
        
        std::cout << "RESULTING SIZE: " << queue_.size() << std::endl;
    }

    void restart(T& n) 
    {
        delete_queue();
        push(n);
    }
    
private:
    own_priority_queue<T, std::vector<T>, CostCompare> queue_;
    std::set<T, FuncCompare>                           membership_;
    unsigned int                                       last_node_copied = 0;
};
#endif
