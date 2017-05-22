/*
 * ibm_apply_runge_kutta.h
 *
 *  Created on: 27 Apr 2017
 *      Author: albert
 */

#ifndef CODE_TEST_IBM_APPLY_RUNGE_KUTTA_H_
#define CODE_TEST_IBM_APPLY_RUNGE_KUTTA_H_

namespace amdados {
namespace app {

template<size_t SizeX, size_t SizeY>
double applyRungeKutta(allscale::utils::grid<double,SizeX,SizeY>& y, double flowu, double flowv, double dt, double ds) {

// TODO
// 1) Add advection part to equations
// 2) compute boundary errors to decide convergence
    using data_type = allscale::utils::grid<double,SizeX,SizeY>;
    flowu = 0; flowv = 0;
    double c = 0.2;
    double diffmult = (dt * c)/pow(ds,2);
    double advmult_u = (flowu * dt)/(2*ds);
    double advmult_v = (flowv * dt)/(2*ds);
    //c = 0.2;
    // solving equation
    //      u_{t+dt}(s) = c * u_{t}(s) / ds^2
    // with linear elements an
    //      u_{t+dt}(s) = c * u_{t}(s) / ds^2

    // the resulting data field
    data_type res;

    // everything but the borders
    for(int i=1; i<(int)SizeX-1; ++i) {      // SizeX and SizeY are number of elements in X and Y
        for(int j=1; j<(int)SizeY-1; ++j) {
            // for some linear element
        	res[{i,j}] =  (y[{i+1,j}] * (diffmult - advmult_u)) +
        			      (y[{i-1,j}] * (diffmult + advmult_u)) +
        			      (y[{i,j+1}] * (diffmult - advmult_v)) +
        			      (y[{i,j-1}] * (diffmult + advmult_v)) +
        			      (y[{i,j  }] * (1 - (4*diffmult)) );
        }
    }

    // upper boarder
    for(int i=0, j=1; j<(int)SizeY-1; ++j) {
        // for some linear element
        res[{i,j}] = y[{i,j}] + c * (y[{i,j}] + y[{i+1,j}] + y[{i,j-1}] + y[{i,j+1}] - 4 * y[{i,j}]) -
            (  (y[{i,j}] - y[{i+1,j}])*flowu/(ds) + (y[{i,j-1}] - y[{i,j+1}])*flowv/(2*ds));;
    }

    // lower boarder
    for(int i=SizeX-1, j=1; j<(int)SizeY-1; ++j) {
        // for some linear element
        res[{i,j}] = y[{i,j}] + c * (y[{i-1,j}] + y[{i,j}] + y[{i,j-1}] + y[{i,j+1}] - 4 * y[{i,j}]) -
            (  (y[{i-1,j}] - y[{i,j}])*flowu/(ds) + (y[{i,j-1}] - y[{i,j+1}])*flowv/(2*ds));
    }

    // left boarder
    for(int i=1,j=0; i<(int)SizeX-1; ++i) {
        // for some linear element
        res[{i,j}] = y[{i,j}] + c * (y[{i-1,j}] + y[{i+1,j}] + y[{i,j}] + y[{i,j+1}] - 4 * y[{i,j}]) -
            (  (y[{i-1,j}] - y[{i+1,j}])*flowu/(2*ds) + (y[{i,j}] - y[{i,j+1}])*flowv/(ds));
    }
    // right boarder
    for(int i=1,j=SizeY-1; i<(int)SizeX-1; ++i) {
        // for some linear element
        res[{i,j}] = y[{i,j}] + c * (y[{i-1,j}] + y[{i+1,j}] + y[{i,j-1}] + y[{i,j}] - 4 * y[{i,j}]) -
            (  (y[{i-1,j}] - y[{i+1,j}])*flowu/(2*ds) + (y[{i,j-1}] - y[{i,j}])*flowv/(ds));
    }

    // and the corners
    res[{      0,      0}] = y[{      0,      0}] + c * (y[{      0,      0}] + y[{      1,      0}] + y[{      0,      0}] + y[{      0,      1}] - 4 * y[{      0,      0}]);
    res[{      0,SizeY-1}] = y[{      0,SizeY-1}] + c * (y[{      0,SizeY-1}] + y[{      1,SizeY-1}] + y[{      0,SizeY-2}] + y[{      0,SizeY-1}] - 4 * y[{      0,SizeY-1}]);
    res[{SizeX-1,      0}] = y[{SizeX-1,      0}] + c * (y[{SizeX-2,      0}] + y[{SizeX-1,      0}] + y[{SizeX-1,      0}] + y[{SizeX-1,      1}] - 4 * y[{SizeX-1,      0}]);
    res[{SizeX-1,SizeY-1}] = y[{SizeX-1,SizeY-1}] + c * (y[{SizeX-2,SizeY-1}] + y[{SizeX-1,SizeY-1}] + y[{SizeX-1,SizeY-2}] + y[{SizeX-1,SizeY-1}] - 4 * y[{SizeX-1,SizeY-1}]);

    // update input
    y = res;

    return 0.0;

    /*
    {
        data_type res;
        bool change = false;
        for(int i=0; i<(int)SizeX; ++i) {
            for(int j=0; j<(int)SizeY; ++j) {
                int sum = 0;
                if (y[{i,j}] > 0) sum = 1;
                if (sum == 0 && i > 0 && y[{i-1,j}] > 0) sum = 1;
                if (sum == 0 && i < SizeX-1 && y[{i+1,j}] > 0) sum = 1;
                if (sum == 0 && j > 0 && y[{i,j-1}] > 0) sum = 1;
                if (sum == 0 && j < SizeY-1 && y[{i,j+1}] > 0) sum = 1;
                res[{i,j}] = sum;
                change = change || res[{j,j}] != y[{i,j}];
            }
        }
        y = res;
        return (change) ? 1.0 : 0.0;
    }
    */

//      // assuming simple, linear stiffness matrix
//      // http://www.mate.tue.nl/~hulsen/cr/slides2.pdf
//
//      // compute k1
    data_type k2;

//      data_type k3;
//      data_type k4;
//
//      computeStep(k1, [&](int i, int j) { return y[{i,j}]; });
//      computeStep(k2, [&](int i, int j) { return y[{i,j}] + delta/2 * k1[{i,j}]; });
//      computeStep(k3, [&](int i, int j) { return y[{i,j}] + delta/2 * k2[{i,j}]; });
//      computeStep(k4, [&](int i, int j) { return y[{i,j}] + delta * k1[{i,j}]; });
//
//      // update input array and compute residuums-norm
//      double res = 0.0;
//      for(int i=0; i<(int)SizeX; ++i) {
//          for(int j=0; j<(int)SizeY; ++j) {
////                y[{i,j}] = y[{i,j}] + delta / 6 * ( k1[{i,j}] + 2*k2[{i,j}] + 2*k3[{i,j}] + k4[{i,j}]);
//              double step = delta / 6 * ( k1[{i,j}] + 2*k2[{i,j}] + 2*k3[{i,j}] + k4[{i,j}]);
//              res += abs(step);
//              y[{i,j}] = y[{i,j}] + step;
//          }
//      }
//
//      return res;
}

template<typename CellConfig>
double applyRungeKutta(Cell<double,CellConfig>& cell, double flowu, double flowv, double delta,double stepsize) {
    if (cell.getActiveLayer() == 0)
    {
    	stepsize = 4;
        return applyRungeKutta(cell.template getLayer<0>(),flowu, flowv, delta,  stepsize);
    }
    if (cell.getActiveLayer() == 1)
    {
    	stepsize = 20;
        return applyRungeKutta(cell.template getLayer<1>(),flowu, flowv, delta, stepsize);
    }
    if (cell.getActiveLayer() == 2)
    {
    	stepsize = 100;
        return applyRungeKutta(cell.template getLayer<2>(),flowu, flowv, delta,stepsize);
    }
//      if (cell.getActiveLayer() == 3)
//          return applyRungeKutta(cell.template getLayer<3>(),delta);
    assert(false && "Not more than 3 layers supported yet!");
    return 0.0;
}


} // namespace app
} // namespace amdados

#endif /* CODE_TEST_IBM_APPLY_RUNGE_KUTTA_H_ */
