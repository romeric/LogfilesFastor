#define FEM_3D_BENCH
#include "helper.h"

#ifdef VECTORISED_CLASSIC_OVERLOADS
    #include "vector_impl.h"
#else
    #include "scalar_impl.h"
#endif


template<typename Ti, typename Tf, size_t ndim, size_t nelem, size_t nnode, size_t nodeperelem, size_t ngauss>
void run(const Tensor<Ti,nelem,nodeperelem> &elements, const Tensor<Tf,nnode,ndim> &points,
        const Tensor<Tf,ndim,nodeperelem,ngauss> &Jm, const Tensor<Tf,ngauss,1> &AllGauss) {

    volatile Tf u1 = random()+2;
    volatile Tf u2 = random()+1;
    volatile Tf ue = random()+0.05;
    volatile Tf e1 = random()+2;
    volatile Tf e2 = random()+1;
    volatile Tf ee = random()+0.05;
    volatile Tf kappa = random()+3;
    auto I = eye<Tf,ndim,ndim,ndim,ndim>();

#ifdef VECTORISED_CLASSIC_OVERLOADS
    const Tf *Jm_data = Jm.data();
#endif

    for (int elem=0; elem<nelem; ++elem) {
        Tensor<Tf,nodeperelem,ndim> LagrangeElemCoords;
        Tensor<Tf,nodeperelem,ndim> EulerElemCoords;
        Tensor<Tf,ndim> D0; D0.random();

        for (int i=0; i<nodeperelem; ++i) {
            for (int j=0; j<ndim; ++j) {
                LagrangeElemCoords(i,j) = points(elements(elem,i),j);
                EulerElemCoords(i,j) = points(elements(elem,i),j);
            }
        }

        for (int g=0; g<ngauss; ++g) {

            Tensor<Tf,ndim,nodeperelem> Jm_g;
            for (int i=0; i<ndim; ++i)
                for (int j=0; j<nodeperelem; ++j)
                    Jm_g(i,j) = Jm(i,j,g);

            // Compute gradient of shape functions
            auto ParentGradientX = sv::matmul(Jm_g,LagrangeElemCoords);
            // Compute material gradient
            auto MaterialGradient = sv::matmul(inverse(ParentGradientX),Jm_g);

            // Compute the deformation gradient tensor
            auto F = sv::matmul(MaterialGradient,EulerElemCoords);
            // Compute H
            auto H = cofactor(F);
            // auto H = sv::cofactor(F);
            // Comput J
            auto J = determinant(F);

            // Compute d=FD0
            auto d = sv::matmul(F,D0);

            // Compute work-conjugates
            Tensor<Tf,ndim,ndim> WF = sv::mul(2.*u1,F) + sv::mul(4*ue*sv::dot(F,F),F) + 2/ee*sv::dot(d,d);
            auto WH = sv::mul(2.*u2,H);
            auto WJ = -(2*u1+4*u2+12*ue)/J+kappa*(J-1) - (1./2./e2/J/J)*sv::dot(d,d);
            auto WD0 = sv::mul(1/e1,D0);
            Tensor<Tf,ndim> Wd = sv::mul(1/e2/J,d) + sv::mul(2/ee*sv::dot(F,F),d)+sv::mul(4/ue/ee/ee*sv::dot(d,d),d);
            // Compute first Piola-Kirchhoff stress tensor
#ifdef OPTIMISED_CROSS
            auto P = sv::add(WF,sv::add(cross(static_cast<Tensor<Tf,ndim,ndim>>(WH),F), sv::mul(WJ,H)));
#else
            auto P = sv::add(WF,sv::add(sv::cross(WH,F), sv::mul(WJ,H)));
#endif
            // Compute electric field
            auto E0 = sv::add(WD0, sv::matmul(transpose(F),Wd));

            // Compute Hessian components
            auto WFF = sv::mul(2.*u1,I) + sv::mul(4*ue*sv::dot(F,F),I)+8*ue*outer(F,F);
            auto WFd = sv::mul(4/ee,d);
            auto WHH = sv::mul(2.*u2,I);
            auto WJJ = 1./J/J*(2*u1+4*u2+12*ue)+kappa + 1/4/e2/J/J/J*sv::dot(d,d);
            auto WJd = sv::mul(-1/e2/J/J,d);
            auto WD0D0 = sv::mul(1/e1,eye2<Tf,ndim,ndim>());
            auto Wdd = sv::mul(1/e2/J,eye2<Tf,ndim,ndim>()) + sv::mul(2/ee*sv::dot(F,F),eye2<Tf,ndim,ndim>()) +
                    sv::mul(4/ue/ee/ee*dot(d,d),eye2<Tf,ndim,ndim>()) +
                    sv::mul(8/ue/ee/ee*dot(d,d),outer(d,d));

            // dump
            unused(P,E0);
            unused(WFF,WFd,WHH,WJJ,WJd,WD0D0,Wdd);
        }
    }
}



FASTOR_NOINLINE void make_run(const std::string &folder) {

    // Obtain quadrature rule, function space and mesh
    std::string efile  = folder+"/meshes/mesh_hand_3d_elements_p"+std::to_string(POLYDEG)+".dat";
    std::string pfile  = folder+"/meshes/mesh_hand_3d_elements_p"+std::to_string(POLYDEG)+".dat";
    std::string gfile0 = folder+"/meshes/p"+std::to_string(POLYDEG)+"_3d_Jm.dat";
    std::string gfile1 = folder+"/meshes/p"+std::to_string(POLYDEG)+"_3d_AllGauss.dat";

    // Load files
    Tensor<size_t,nelem,nodeperelem> elements = loadtxt<size_t,nelem,nodeperelem>(efile);
    Tensor<real,nnode,ndim> points = loadtxt<real,nnode,ndim>(pfile);

    Tensor<real,ndim*nodeperelem*ngauss,1> Jm_temp = loadtxt<real,ndim*nodeperelem*ngauss,1>(gfile0);
    Tensor<real,ndim,nodeperelem,ngauss> Jm;
    std::copy(Jm_temp.data(),Jm_temp.data()+Jm_temp.Size,Jm.data()); 
    Tensor<real,ngauss,1> AllGauss = loadtxt<real,ngauss,1>(gfile1);

    // warm-up call
    run(elements,points,Jm,AllGauss);

    // run benchmark
    double time;
    std::tie(time,std::ignore) = rtimeit(static_cast<void (*)(const Tensor<size_t,nelem,nodeperelem> &, const Tensor<real,nnode,ndim> &,
        const Tensor<real,ndim,nodeperelem,ngauss> &, const Tensor<real,ngauss,1> &)>(&run),elements,points,Jm,AllGauss);
    print(time);
}



int main(int argc, char *argv[]) {

    FASTOR_ASSERT(argc!=1, BOLD(FBLU("Number of arguments to the main function should be 1 (i.e. name of the FEM folder)")));
    std::string folder = argv[1];


    //--------------------------------------------------------------------------------//
    const rlim_t stacksize = 1024*1024*1024;
    struct rlimit rl;
    int result;
    result = getrlimit(RLIMIT_STACK, &rl);
    if (result==0) {
        if (rl.rlim_cur < stacksize) {
            rl.rlim_cur = stacksize;
            result = setrlimit(RLIMIT_STACK,&rl);
            if (result !=0) {
                FASTOR_ASSERT(result !=0, "CHANGING STACK SIZE FAILED");
            }
        }
    }
   //--------------------------------------------------------------------------------//

    make_run(folder);

    return 0;
}