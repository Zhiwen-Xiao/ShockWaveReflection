/*
 * @version: v0.1
 * @Descripttion: 给出一些基本的激波关系式和求解函数，利用了自编的Solver求解器
 * @Author: XiaoZhiwen
 * @Date: 2019-11-09 17:13:12
 */
#ifndef __SHOCKWAVE__HPP__
#define __SHOCKWAVE__HPP__

#include "inc.hpp"
#include "Solver.hpp"

typedef long double REAL;

#define Pi 3.141592653589793115997963468544185161590576171875

/**
 * @ description 描述: 已知来流马赫数Ma、激波角 beta 求气流偏转角 theta
 * @ param 参数: 
 *      Ma      :来流马赫数
 *      beta    :激波角（弧度制）
 *      gamma   :比热比，默认1.4
 * @ return 返回值: 
 *      TanTheta:气流偏转角
 */  
REAL Theta(REAL beta, REAL Ma, REAL gamma = 1.4)
{  
    REAL tantheta = 2 * (pow(Ma,2) * pow(sin(beta),2) - 1) / (tan(beta) * (pow(Ma,2) * (gamma + cos(2 * beta)) + 2));
    return atan(tantheta);
}

/**
 * @ description 描述: 已知来流马赫数Ma、激波角 beta 求激波前后气压比P2/P1
 * @ param 参数: 
 *      Ma      :来流马赫数
 *      beta    :激波角（弧度制）
 *      gamma   :比热比，默认1.4
 * @ return 返回值: 
 *      P2/P1
 */
REAL P_Ratio(REAL Ma, REAL beta = Pi/2, REAL gamma = 1.4)
{
    return 1 + (2 * gamma) / (gamma + 1) * (pow(Ma,2) * pow(sin(beta),2) - 1);
}

/**
 * @ description 描述: 已知来流马赫数Ma、激波角 beta 求激波前后密度比rho2/rho1
 * @ param 参数: 
 *      Ma      :来流马赫数
 *      beta    :激波角（弧度制）
 *      gamma   :比热比，默认1.4
 * @ return 返回值: 
 *      rho2/rho1
 */
REAL Rho_Ratio(REAL Ma, REAL beta = Pi/2, REAL gamma = 1.4)
{
    REAL temp = pow(Ma,2) * pow(sin(beta),2);
    return ((gamma + 1) * temp) / (2 + (gamma - 1) * temp);
}

/**
 * @ description 描述: 已知来流马赫数Ma、激波角 beta 求激波前后温度比T2/T1
 * @ param 参数: 
 *      Ma      :来流马赫数
 *      beta    :激波角（弧度制）
 *      gamma   :比热比，默认1.4
 * @ return 返回值: 
 *      T2/T1
 */
REAL T_Ratio(REAL Ma, REAL beta = Pi/2, REAL gamma = 1.4)
{
    REAL temp = pow(Ma,2) * pow(sin(beta),2);
    return (2 * gamma * temp - (gamma - 1)) * ((gamma - 1) * temp + 2) / (pow((gamma + 1),2) * temp);
}

/**
 * @ description 描述: 已知来流马赫数Ma、激波角 beta 求激波后马赫数
 * @ param 参数: 
 *      Ma      :来流马赫数
 *      beta    :激波角（弧度制）
 *      gamma   :比热比，默认1.4
 * @ return 返回值: 
 *      Ma2
 */
REAL Ma2(REAL Ma, REAL beta = Pi/2, REAL gamma = 1.4)
{
    REAL temp1 = pow(Ma,2) * pow(sin(beta),2);
    REAL temp2 = pow(Ma,2) * pow(cos(beta),2);
    REAL Ma2 = (pow(Ma,2) + 2 / (gamma - 1)) / ((2 * gamma) / (gamma - 1) * temp1 - 1) + temp2 / ((gamma - 1) / 2 * temp1 + 1);

    return sqrt(Ma2);
}

/**
 * @ description 描述: 已知来流马赫数Ma，求最大气流偏转角对应激波角 betaMax .
 * @ param 参数: 
 *      Ma      :来流马赫数
 *      gamma   :比热比，默认1.4
 * @ return 返回值: 
 *      betaMax :最大气流偏转角对应激波角
 */  
REAL BetaMax(REAL Ma, REAL gamma = 1.4)
{  
    REAL temp = (gamma + 1) * (1 + (gamma - 1) / 2 * pow(Ma,2) + (gamma + 1) / 16 * pow(Ma,4));
    REAL sinbeta2 = 1 / (gamma * pow(Ma,2)) * ((gamma + 1) / 4 * pow(Ma,2) - 1 + sqrt(temp));
    return asin(sqrt(sinbeta2));
}

/**
 * @ description 描述: 已知来流马赫数Ma，求斜激波脱体角 thetaMax .
 * @ param 参数: 
 *      Ma      :来流马赫数
 *      gamma   :比热比，默认1.4
 * @ return 返回值: 
 *      thetaMax :最大气流偏转角对应激波角
 */  
REAL ThetaMax(REAL Ma, REAL gamma = 1.4)
{  
    REAL betaMax = BetaMax(Ma, gamma);
    return Theta(betaMax, Ma, gamma);
}

/**
 * @ description 描述: 已知来流马赫数Ma、气流偏转角 theta，使用拟牛顿法求解激波角 beta，
 * @ param 参数: 
 *      Ma      :来流马赫数
 *      beta    :激波角（弧度制）
 *      option  :强/弱解选择：0: 弱解； 1: 强解
 *      gamma   :比热比，默认1.4
 * @ return 返回值: 
 *      TanTheta:气流偏转角
 */  
REAL Beta(REAL theta, REAL Ma, int option = 0, REAL gamma = 1.4)
{  
    REAL thetaMax = ThetaMax(Ma,gamma);
    if (theta >= thetaMax)
    {
        return Pi/2;
    } else{
        using namespace std::placeholders;
        auto thetaFun = std::bind (Theta, _1, Ma, gamma);
        REAL beta;
        REAL betaMax = BetaMax(Ma,gamma);
        if (option == 0)
        {
            beta = QuasiNewtonSolver(thetaFun,0.1 * betaMax,0.9 * betaMax,theta);
        }else
        {
            beta = QuasiNewtonSolver(thetaFun,betaMax + 0.1 * (Pi/2 - betaMax),betaMax + 0.9 * (Pi/2 - betaMax),theta);
        }
        
        
        return beta;
    }
}

/**
 * @ description 描述: 激波反射中，已知来流马赫数Ma、第一个楔角（气流转折角）theta，求2区脱体角 .
 * @ param 参数: 
 *      theta   :气流折转角
 *      Ma      :来流马赫数
 *      gamma   :比热比，默认1.4
 * @ return 返回值: 
 *      thetaMax2 :2区脱体角
 */  
REAL ThetaMax2(REAL theta, REAL Ma, REAL gamma = 1.4)
{  
    REAL beta1 = Beta(theta, Ma, 0, gamma);
    REAL ma2 = Ma2(Ma, beta1, gamma);
    REAL thetaMax2 = ThetaMax(ma2, gamma);
    return thetaMax2;
}

/**
 * @ description 描述: 已知来流马赫数Ma，使用SteffensenSolver求解马赫反射脱体条件角 MachConditionTheta .
 * @ param 参数: 
 *      Ma      :来流马赫数
 *      gamma   :比热比，默认1.4
 * @ return 返回值: 
 *      machConditiontheta :马赫反射脱体条件角
 */  
REAL MachConditionTheta(REAL Ma, REAL gamma = 1.4)
{  
    using namespace std::placeholders;
    auto thetaFun = std::bind (ThetaMax2, _1, Ma, gamma);
    REAL machConditionTheta = SteffensenSolver(thetaFun,0.01);
    return machConditionTheta;
}

/**
 * @ description 描述: 激波反射中，已知来流马赫数Ma、第一个楔角（气流转折角）theta，求正规反射下2区亚比P2/P0 .
 * @ param 参数: 
 *      theta   :气流折转角
 *      Ma      :来流马赫数
 *      gamma   :比热比，默认1.4
 * @ return 返回值: 
 *      P2/P0   :2区
 */  
REAL P2_Ratio_N(REAL theta, REAL Ma, REAL gamma = 1.4)
{  
    REAL beta1 = Beta(theta, Ma, 0, gamma);
    REAL p1 = P_Ratio(Ma, beta1, gamma);
    REAL ma1 = Ma2(Ma, beta1, gamma);
    REAL beta2 = Beta(theta, ma1, 0, gamma);
    REAL p2 = P_Ratio(ma1, beta2, gamma);
    return p2 * p1;
}

/**
 * @ description 描述: 已知来流马赫数Ma，使用QuasiNewtonSolver求解马赫反射von-Neumann压力平衡条件角 VonNeumannConditionTheta .
 * @ param 参数: 
 *      Ma      :来流马赫数
 *      gamma   :比热比，默认1.4
 * @ return 返回值: 
 *      machConditiontheta :马赫反射脱体条件角
 */  
REAL VonNeumannConditionTheta(REAL Ma, REAL gamma = 1.4)
{  
    REAL pn = P_Ratio(Ma);
    REAL theta0 = MachConditionTheta(Ma, gamma);

    using namespace std::placeholders;
    auto pFun = std::bind (P2_Ratio_N, _1, Ma, gamma);
    REAL vonNeumannConditionTheta = QuasiNewtonSolver(pFun,theta0,theta0/2,pn);
    return vonNeumannConditionTheta;
}

/**
 * @ description 描述: 激波反射中，已知来流马赫数Ma、第一个楔角（气流转折角）theta，滑移角 thetaS, 求马赫反射下2区亚比P2/P0与3区亚比P3/P0之差 .
 * @ param 参数: 
 *      thetaS  :滑移线与水平线夹角
 *      theta   :气流折转角
 *      Ma      :来流马赫数
 *      gamma   :比热比，默认1.4
 * @ return 返回值: 
 *      P2/P0-P3/P0   :2区
 */  
REAL DeltaP23_Ratio_M(REAL thetaS, REAL theta, REAL Ma, REAL gamma = 1.4)
{  
    REAL beta1 = Beta(theta, Ma, 0, gamma);
    REAL p1 = P_Ratio(Ma, beta1, gamma);
    REAL ma1 = Ma2(Ma, beta1, gamma);
    REAL beta2 = Beta(theta-thetaS, ma1, 0, gamma);
    REAL p2 = P_Ratio(ma1, beta2, gamma);
    p2 = p1 * p2;
    REAL beta3 = Beta(thetaS,Ma,1,gamma);
    REAL p3 = P_Ratio(Ma, beta3, gamma);

    return (p2 - p3) * 100000;
}

/**
 * @ description 描述: 已知来流马赫数Ma、第一个楔角（气流转折角）theta，使用QuasiNewtonSolver求解马赫反射滑移角 thetaS.
 * @ param 参数: 
 *      theta   :气流折转角
 *      Ma      :来流马赫数
 *      gamma   :比热比，默认1.4
 * @ return 返回值: 
 *      machConditiontheta :马赫反射脱体条件角
 */  
REAL ThetaS_M(REAL theta, REAL Ma, REAL gamma = 1.4)
{  
    REAL thetaMax2 = ThetaMax2(theta, Ma, gamma);
    REAL deltaTheta = abs(thetaMax2 - theta);
    using namespace std::placeholders;
    auto deltaPFun = std::bind (DeltaP23_Ratio_M, _1, theta, Ma, gamma);
    REAL thetaS = QuasiNewtonSolver(deltaPFun,0.1 * deltaTheta, 2 * deltaTheta);
    return thetaS;
}


#endif