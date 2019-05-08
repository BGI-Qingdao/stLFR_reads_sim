#ifndef __APPCOMMON_POOL_H__
#define __APPCOMMON_POOL_H__

#include <malloc.h>

namespace BGIQD {
    namespace stLFRSim {

        template< class U >
            struct Pool
            {
                public:
                    typedef U Unit;
                    Pool() 
                        : unit_size(sizeof(Unit)) 
                          , data(NULL) 
                          , curr_size(1024)
                          , curr_top(1024)
                    {
                        data = static_cast<Unit*>( malloc( curr_size * unit_size ) );
                        assert(data);
                    }

                    Unit * Push() 
                    {
                        if( curr_top < curr_size )
                            return &(data[curr_top++]);
                        curr_size *= 2 ;
                        data =static_cast<Unit*>( realloc(data,curr_size * unit_size) );
                        assert(data);
                        assert(curr_top < curr_size );
                        return &(data[curr_top++]);
                    }
                    Unit * Top() const
                    {
                        assert(curr_top > 0);
                        assert(curr_top <= curr_size );
                        return &(data[curr_top-1]);
                    }
                    void Pop()
                    {
                        assert(curr_top > 0);
                        curr_top -- ;
                    }
                    int Size() const { return curr_top ; }
                    void Clear()
                    {
                        curr_top = 0 ;
                    }
                private:
                    const size_t unit_size ;
                    Unit * data;
                    int curr_size ;
                    int curr_top ;
            };
    }
}

#endif // __APPCOMMON_POOL_H__
