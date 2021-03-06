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
    volatile Tf kappa = random()+3;
    auto I = eye<Tf,ndim,ndim,ndim,ndim>();

#ifdef VECTORISED_CLASSIC_OVERLOADS
    const Tf *Jm_data = Jm.data();
#endif

    for (int elem=0; elem<nelem; ++elem) {
        Tensor<Tf,nodeperelem,ndim> LagrangeElemCoords;
        Tensor<Tf,nodeperelem,ndim> EulerElemCoords;

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

            // Compute work-conjugates
            auto WF = sv::mul(2.*u1,F);
            auto WH = sv::mul(2.*u2,H);
            auto WJ = -2.*(u1+2*u2)/J+kappa*(J-1);

            // Compute first Piola-Kirchhoff stress tensor
#ifdef OPTIMISED_CROSS
            auto P = sv::add(WF,sv::add(cross<PlaneStrain>(static_cast<Tensor<Tf,ndim,ndim>>(WH),F), sv::mul(WJ,H)));
#else
            auto P = sv::add(WF,sv::add(sv::cross(static_cast<Tensor<Tf,ndim,ndim>>(WH),F), sv::mul(WJ,H)));
#endif
            // Compute Hessian components
            auto WFF = sv::mul(2.*u1,I);
            auto WHH = sv::mul(2.*u2,I);
            auto WJJ = 2./J/J*(u1+2*u2)+kappa;

            // dump
            unused(P);
            unused(WFF,WHH,WJJ);
        }
    }
}

int main(int argc, char *argv[]) {

    FASTOR_ASSERT(argc!=1, BOLD(FBLU("Number of arguments to the main function should be 1 (i.e. name of the FEM folder)")));
    std::string folder = argv[1];

    // Obtain quadrature rule, function space and mesh
    std::string efile  = folder+"/meshes/mesh_2d_elements_p"+std::to_string(POLYDEG)+".dat";
    std::string pfile  = folder+"/meshes/mesh_2d_points_p"+std::to_string(POLYDEG)+".dat";
    std::string gfile0 = folder+"/meshes/p"+std::to_string(POLYDEG)+"_2d_Jm.dat";
    std::string gfile1 = folder+"/meshes/p"+std::to_string(POLYDEG)+"_2d_AllGauss.dat";

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
    return 0;
}



