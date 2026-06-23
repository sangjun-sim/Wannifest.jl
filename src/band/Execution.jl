module Execution

using LinearAlgebra: BLAS

export with_blas_threads

function with_blas_threads(f::Function, num_threads::Integer)
    num_threads > 0 || error("BLAS thread count must be positive")
    old_threads = BLAS.get_num_threads()
    old_threads == num_threads && return f()
    BLAS.set_num_threads(num_threads)
    try
        return f()
    finally
        BLAS.set_num_threads(old_threads)
    end
end

end
