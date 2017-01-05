#define FEM_3D_BENCH
#include "helper.h"


template<typename Ti, typename Tf, size_t ndim, size_t nelem, size_t nnode, size_t nodeperelem, size_t ngauss>
void run(const Tensor<Ti,nelem,nodeperelem> &elements, const Tensor<Tf,nnode,ndim> &points,
        const Tensor<Tf,ndim,nodeperelem,ngauss> &Jm, const Tensor<Tf,ngauss,1> &AllGauss) {

    volatile Tf u1 = random()+2;
    volatile Tf u2 = random()+1;
    volatile Tf e1 = random()+2;
    volatile Tf e2 = random()+1;
    volatile Tf kappa = random()+3;
    auto I = eye<Tf,ndim,ndim,ndim,ndim>();

    const Tf *Jm_data = Jm.data();

    for (int elem=0; elem<nelem; ++elem) {
        Tensor<Tf,nodeperelem,ndim> LagrangeElemCoords;
        Tensor<Tf,nodeperelem,ndim> EulerElemCoords;
        Tensor<Tf,ndim> D0; D0.random();

        for (int i=0; i<nodeperelem; ++i) {
            for (int j=0; j<ndim; ++j) {
                LagrangeElemCoords(i,j) = points(elements(elem,i),j);
                EulerElemCoords(i,j) = points(elements(elem,i),j) + random();
            }
        }

        for (int g=0; g<ngauss; ++g) {
            // Get Gauss point Jm
            Tensor<Tf,ndim,nodeperelem> Jm_g;
            for (int i=0; i<ndim; ++i)
                for (int j=0; j<nodeperelem; ++j)
                    Jm_g(i,j) = Jm(i,j,g);

            auto ParentGradientX = lmatmul(Jm_g,LagrangeElemCoords);
            // Compute inverse - needs to be lvalue otherwise the smart expression is invalid
            auto invParentGradientX = linverse(ParentGradientX);
            // Compute material gradient
            auto MaterialGradient = lmatmul(invParentGradientX,Jm_g);

            // Compute the deformation gradient tensor
            Tensor<Tf,ndim,ndim> F = lmatmul(MaterialGradient,EulerElemCoords);

            // Compute H
            auto H = cofactor(F);
            // Comput J
            auto J = determinant(F);

            // Compute d=FD0
            auto d = matmul(F,D0);

            // Compute work-conjugates
            auto WF = 2.*u1*F;
            auto WH = 2.*u2*H;
            auto WJ = -2.*(u1+2*u2)/J+kappa*(J-1) - (1./2./e2/J/J)*dot(d,d);
            auto WD0 = 1/e1*D0;
            auto Wd = 1/e2/J*d;
            // Compute first Piola-Kirchhoff stress tensor
            auto P = WF + cross(static_cast<Tensor<Tf,ndim,ndim>>(WH),F) + WJ*H;
            // Compute electric field
            auto E0 = WD0 + matmul(transpose(F),static_cast<Tensor<Tf,ndim>>(Wd));

            // Compute Hessian components
            auto WFF = 2.*u1*I;
            auto WHH = 2.*u2*I;
            auto WJJ = 2./J/J*(u1+2*u2)+kappa + 1/4/e2/J/J/J*dot(d,d);
            auto WJd = -1/e2/J/J*d;
            auto WD0D0 = 1/e1*eye2<Tf,ndim,ndim>();
            auto Wdd = 1/e2/J*eye2<Tf,ndim,ndim>();

            // dump
            unused(P,E0);
            unused(WFF,WHH,WJJ,WJd,WD0D0,Wdd);

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