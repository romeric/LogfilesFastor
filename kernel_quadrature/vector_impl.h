#ifndef VECTOR_IMPL_H
#define VECTOR_IMPL_H

#include "helper.h"

// sv: vector variant also uses the same namespace
namespace sv {

template<typename T, size_t M, size_t N, 
    typename std::enable_if<M==3 && N==3,bool>::type=0>
inline Tensor<T,M,N> cofactor(const Tensor<T,M,N> &a) {
    constexpr T levi_civita[27] = { 0.,  0.,  0.,  0.,  0., -1.,  0.,  1.,  0.,  0.,  0.,  1.,  0.,  0.,
                                0., -1.,  0.,  0.,  0., -1.,  0.,  1.,  0.,  0.,  0.,  0.,  0.};
    Tensor<T,M,N> out;
    constexpr size_t size = N;
    for (size_t i=0; i<N; ++i)
        for (size_t j=0; j<N; ++j)
            for (size_t k=0; k<N; ++k)
                for (size_t I=0; I<N; ++I)
                    for (size_t J=0; J<N; ++J)
                        for (size_t K=0; K<N; ++K)
                            out(i,I) += 0.5*levi_civita[i*size*size+j*size+k]*levi_civita[I*size*size+J*size+K]*a(j,J)*a(k,K);
    return out;
}

template<typename T, size_t M, size_t N, 
    typename std::enable_if<M==2 && N==2,bool>::type=0>
inline Tensor<T,M,N> cofactor(const Tensor<T,M,N> &a) {
    constexpr T levi_civita[27] = { 0.,  0.,  0.,  0.,  0., -1.,  0.,  1.,  0.,  0.,  0.,  1.,  0.,  0.,
                                0., -1.,  0.,  0.,  0., -1.,  0.,  1.,  0.,  0.,  0.,  0.,  0.};

    Tensor<T,3,3> a3d; a3d.zeros();

    for (size_t i=0; i<2; ++i) {
        for (size_t j=0; j<2; ++j) {
            a3d(i,j) = a(i,j);
        }
    }
    a3d(2,2) = 1;

    Tensor<T,M+1,N+1> out;
    constexpr size_t size = N+1;
    for (size_t i=0; i<N; ++i)
        for (size_t j=0; j<N; ++j)
            for (size_t k=0; k<N; ++k)
                for (size_t I=0; I<N; ++I)
                    for (size_t J=0; J<N; ++J)
                        for (size_t K=0; K<N; ++K)
                            out(i,I) += 0.5*levi_civita[i*size*size+j*size+k]*levi_civita[I*size*size+J*size+K]*a3d(j,J)*a3d(k,K);

    Tensor<T,M,N> out2d;
    for (size_t i=0; i<2; ++i) {
        for (size_t j=0; j<2; ++j) {
            out2d(i,j) = out(i,j);
        }
    }

    return out2d;
}


template<typename T, size_t M, size_t N, 
    typename std::enable_if<M==3 && N==3,bool>::type=0>
inline Tensor<T,M,N> cross(const Tensor<T,M,N> &a, const Tensor<T,M,N> &b) {
    constexpr T levi_civita[27] = { 0.,  0.,  0.,  0.,  0., -1.,  0.,  1.,  0.,  0.,  0.,  1.,  0.,  0.,
                                0., -1.,  0.,  0.,  0., -1.,  0.,  1.,  0.,  0.,  0.,  0.,  0.};
    Tensor<T,M,N> out;
    constexpr size_t size = N;
    for (size_t i=0; i<N; ++i)
        for (size_t j=0; j<N; ++j)
            for (size_t k=0; k<N; ++k)
                for (size_t I=0; I<N; ++I)
                    for (size_t J=0; J<N; ++J)
                        for (size_t K=0; K<N; ++K)
                            out(i,I) += levi_civita[i*size*size+j*size+k]*levi_civita[I*size*size+J*size+K]*a(j,J)*b(k,K);
    return out;
}

template<typename T, size_t M, size_t N, 
    typename std::enable_if<M==2 && N==2,bool>::type=0>
inline Tensor<T,M,N> cross(const Tensor<T,M,N> &a, const Tensor<T,M,N> &b) {
    constexpr T levi_civita[27] = { 0.,  0.,  0.,  0.,  0., -1.,  0.,  1.,  0.,  0.,  0.,  1.,  0.,  0.,
                                0., -1.,  0.,  0.,  0., -1.,  0.,  1.,  0.,  0.,  0.,  0.,  0.};

    Tensor<T,3,3> a3d; a3d.zeros();
    Tensor<T,3,3> b3d; b3d.zeros();

    for (size_t i=0; i<2; ++i) {
        for (size_t j=0; j<2; ++j) {
            a3d(i,j) = a(i,j);
            b3d(i,j) = b(i,j);
        }
    }
    a3d(2,2) = 1;
    b3d(2,2) = 1;

    Tensor<T,M+1,N+1> out;
    constexpr size_t size = N+1;
    for (size_t i=0; i<N; ++i)
        for (size_t j=0; j<N; ++j)
            for (size_t k=0; k<N; ++k)
                for (size_t I=0; I<N; ++I)
                    for (size_t J=0; J<N; ++J)
                        for (size_t K=0; K<N; ++K)
                            out(i,I) += levi_civita[i*size*size+j*size+k]*levi_civita[I*size*size+J*size+K]*a3d(j,J)*b3d(k,K);

    Tensor<T,M,N> out2d;
    for (size_t i=0; i<2; ++i) {
        for (size_t j=0; j<2; ++j) {
            out2d(i,j) = out(i,j);
        }
    }

    return out2d;
}


template<typename T, size_t M, size_t N, size_t K>
inline Tensor<T,M,N> matmul(const Tensor<T,M,K> &a, const Tensor<T,K,N> &b) {
    Tensor<T,M,N> out;
    unused(a,b);
    return Fastor::matmul(a,b);
}

template<typename T, size_t M, size_t K>
inline Tensor<T,M> matmul(const Tensor<T,M,K> &a, const Tensor<T,K> &b) {
    Tensor<T,M> out;
    out = Fastor::matmul(a,b);
    unused(out);
    return out;
}

template<typename T, size_t ... Rest>
inline T dot(const Tensor<T,Rest...> &a, const Tensor<T,Rest...> &b) {
    T out = Fastor::dot(a,b);
    unused(a,b,out);
    return out;
}


// add
template<typename T, size_t ... Rest>
inline Tensor<T,Rest...> add(const Tensor<T,Rest...> &a, const Tensor<T,Rest...> &b) {
    Tensor<T,Rest...> out;
    const T *a_data = a.data();
    const T *b_data = b.data();
    T *out_data = out.data();
    using V = SIMDVector<T>;
    V vec;
    for (auto i=0; i<a.Size; i+=V::Size) {
        vec = V(&a_data[i])+V(&b_data[i]);
        vec.store(&out_data[i]);
    }
    unused(a,b,out);
    return out;
}

template<typename T, size_t ... Rest>
inline Tensor<T,Rest...> add(T num, const Tensor<T,Rest...> &b) {
    Tensor<T,Rest...> out;
    const T *b_data = b.data();
    T *out_data = out.data();
    using V = SIMDVector<T>;
    V veca = num;
    V vec;
    for (FASTOR_VINDEX i=0; i<b.Size; i+=V::Size) {
        vec = veca+V(&b_data[i]);
        vec.store(&out_data[i]);
    }
    return out;
}

template<typename T, size_t ... Rest>
inline Tensor<T,Rest...> add(const Tensor<T,Rest...> &b, T num) {
    Tensor<T,Rest...> out;
    const T *b_data = b.data();
    T *out_data = out.data();
    using V = SIMDVector<T>;
    V veca = num;
    V vec;
    for (FASTOR_VINDEX i=0; i<b.Size; i+=V::Size) {
        vec = veca+V(&b_data[i]);
        vec.store(&out_data[i]);
    }
    return out;
}

// sub
template<typename T, size_t ... Rest>
inline Tensor<T,Rest...> sub(const Tensor<T,Rest...> &a, const Tensor<T,Rest...> &b) {
    Tensor<T,Rest...> out;
    const T *a_data = a.data();
    const T *b_data = b.data();
    T *out_data = out.data();
    using V = SIMDVector<T>;
    V vec;
    for (auto i=0; i<a.Size; i+=V::Size) {
        vec = V(&a_data[i])-V(&b_data[i]);
        vec.store(&out_data[i]);
    }
    return out;
}

template<typename T, size_t ... Rest>
inline Tensor<T,Rest...> sub(T num, const Tensor<T,Rest...> &b) {
    Tensor<T,Rest...> out;
    const T *b_data = b.data();
    T *out_data = out.data();
    using V = SIMDVector<T>;
    V veca = num;
    V vec;
    for (auto i=0; i<b.Size; i+=V::Size) {
        vec = veca-V(&b_data[i]);
        vec.store(&out_data[i]);
    }
    return out;
}

template<typename T, size_t ... Rest>
inline Tensor<T,Rest...> sub(const Tensor<T,Rest...> &b, T num) {
    Tensor<T,Rest...> out;
    const T *b_data = b.data();
    T *out_data = out.data();
    using V = SIMDVector<T>;
    V veca = num;
    V vec;
    for (auto i=0; i<b.Size; i+=V::Size) {
        vec = V(&b_data[i]) - veca;
        vec.store(&out_data[i]);
    }
    return out;
}

// mul
template<typename T, size_t ... Rest>
inline Tensor<T,Rest...> mul(const Tensor<T,Rest...> &a, const Tensor<T,Rest...> &b) {
    Tensor<T,Rest...> out;
    const T *a_data = a.data();
    const T *b_data = b.data();
    T *out_data = out.data();
    using V = SIMDVector<T>;
    V vec;
    for (auto i=0; i<a.Size; i+=V::Size) {
        vec = V(&a_data[i])*V(&b_data[i]);
        vec.store(&out_data[i]);
    }
    return out;
}

template<typename T, size_t ... Rest>
inline Tensor<T,Rest...> mul(T num, const Tensor<T,Rest...> &b) {
    Tensor<T,Rest...> out;
    const T *b_data = b.data();
    T *out_data = out.data();
    using V = SIMDVector<T>;
    V veca = num;
    V vec;
    for (FASTOR_VINDEX i=0; i<b.Size; i+=V::Size) {
        vec = veca*V(&b_data[i]);
        vec.store(&out_data[i]);
    }
    return out;
}

template<typename T, size_t ... Rest>
inline Tensor<T,Rest...> mul(const Tensor<T,Rest...> &b, T num) {
    Tensor<T,Rest...> out;
    const T *b_data = b.data();
    T *out_data = out.data();
    using V = SIMDVector<T>;
    V veca = num;
    V vec;
    for (FASTOR_VINDEX i=0; i<b.Size; i+=V::Size) {
        vec = V(&b_data[i])*veca;
        vec.store(&out_data[i]);
    }
    return out;
}

// div
template<typename T, size_t ... Rest>
inline Tensor<T,Rest...> div(const Tensor<T,Rest...> &a, const Tensor<T,Rest...> &b) {
    Tensor<T,Rest...> out;
    const T *a_data = a.data();
    const T *b_data = b.data();
    T *out_data = out.data();
    using V = SIMDVector<T>;
    V vec;
    for (auto i=0; i<a.Size; i+=V::Size) {
        vec = V(&a_data[i])/V(&b_data[i]);
        vec.store(&out_data[i]);
    }
    return out;
}

template<typename T, size_t ... Rest>
inline Tensor<T,Rest...> div(T num, const Tensor<T,Rest...> &b) {
    Tensor<T,Rest...> out;
    const T *b_data = b.data();
    T *out_data = out.data();
    using V = SIMDVector<T>;
    V veca = num;
    V vec;
    for (auto i=0; i<b.Size; i+=V::Size) {
        vec = veca/V(&b_data[i]);
        vec.store(&out_data[i]);
    }
    return out;
}

template<typename T, size_t ... Rest>
inline Tensor<T,Rest...> div(const Tensor<T,Rest...> &b, T num) {
    Tensor<T,Rest...> out;
    const T *b_data = b.data();
    T *out_data = out.data();
    using V = SIMDVector<T>;
    V veca = num;
    V vec;
    for (auto i=0; i<b.Size; i+=V::Size) {
        vec = V(&b_data[i])/veca;
        vec.store(&out_data[i]);
    }
    return out;
}


}


#endif // VECTOR_IMPL_H