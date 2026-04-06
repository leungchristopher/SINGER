//
//  MemoryPool.hpp
//  ChrisSINGER
//
//  Fixed-block pool allocator for Interval objects.
//

#ifndef MemoryPool_hpp
#define MemoryPool_hpp

#include <cstddef>
#include <vector>
#include <cassert>

template <typename T, size_t BlockSize = 4096>
class MemoryPool {

public:

    MemoryPool() = default;

    ~MemoryPool() {
        for (char *block : blocks_) {
            ::operator delete(block);
        }
    }

    MemoryPool(const MemoryPool &) = delete;
    MemoryPool &operator=(const MemoryPool &) = delete;

    T *allocate() {
        if (free_list_) {
            T *ptr = free_list_;
            free_list_ = *reinterpret_cast<T **>(ptr);
            return ptr;
        }
        if (current_offset_ >= BlockSize) {
            allocate_block();
        }
        T *ptr = reinterpret_cast<T *>(current_block_ + current_offset_ * slot_size_);
        ++current_offset_;
        return ptr;
    }

    void deallocate(T *ptr) {
        if (ptr == nullptr) return;
        ptr->~T();
        *reinterpret_cast<T **>(ptr) = free_list_;
        free_list_ = ptr;
    }

    template <typename... Args>
    T *construct(Args &&...args) {
        T *ptr = allocate();
        new (ptr) T(std::forward<Args>(args)...);
        return ptr;
    }

    void reset() {
        free_list_ = nullptr;
        current_offset_ = BlockSize; // force new block on next allocate
    }

private:

    static constexpr size_t slot_size_ = sizeof(T) < sizeof(T *) ? sizeof(T *) : sizeof(T);

    std::vector<char *> blocks_;
    char *current_block_ = nullptr;
    size_t current_offset_ = BlockSize;
    T *free_list_ = nullptr;

    void allocate_block() {
        char *block = static_cast<char *>(::operator new(BlockSize * slot_size_));
        blocks_.push_back(block);
        current_block_ = block;
        current_offset_ = 0;
    }
};

#endif /* MemoryPool_hpp */
