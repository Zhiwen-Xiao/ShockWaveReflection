/*
 * @version: v0.1
 * @Descripttion: 给出两种格式的非线性方程求解器：
 *      SteffensenSolver    可以用来求解迭代格式：x == f(x) 的方程，使用Steffensen加速算法，二阶收敛
 *      QuasiNewtonSolver   可以用来求解等式格式：f(x) == a 的方程，使用拟牛顿算法，超线性收敛
 * @Author: XiaoZhiwen
 * @Date: 2019-11-09 17:20:21
 */
#ifndef __SOLVER__HPP__
#define __SOLVER__HPP__

#include "inc.hpp"

typedef long double REAL;

/**
 * @ description 描述: Steffensen 迭代公式，至少二阶收敛，注意这里求解的是迭代方程x=format(x)
 * @ param 参数: 
 *      format  : 迭代公式函数名
 *      x       : 变量
 * @ return 返回值: 
 *      迭代单次
 */
REAL SteffensenFormat(std::function<REAL(REAL)> format,REAL x)
{
    REAL x1 = format(x);
    REAL x2 = format(x1);
    return x-pow(x1-x,(REAL)2)/(x2-2*x1 +x);
}

/**
 * @ description 描述: Steffensen 迭代法求解器，至少二阶收敛，注意这里求解的是迭代方程 x=format(x)
 * @ param 参数: 
 *      format      : 迭代公式函数名
 *      initValue   : 迭代初始值
 *      errMax      : 最大误差限
 *      MaxLoop     : 最大迭代次数
 * @ return 返回值: 
 *      求解结果
 */
REAL SteffensenSolver(std::function<REAL(REAL)> format, REAL initValue, REAL errMax = 1E-12, int MaxLoop = 20)
{
    REAL xk = initValue;
    REAL xkk = SteffensenFormat(format,xk);
    REAL err = abs(xkk-xk);
    int loop = 1;
    while (loop < MaxLoop && err > errMax)
    {
        xk = xkk;
        xkk = SteffensenFormat(format,xk);
        err = abs(xkk-xk);
        loop ++;
    }
    return xkk;
}

/**
 * @ description 描述: 拟牛顿法求解器，超线性收敛，直接求解方程 format(x) == a
 * @ param 参数:  
     *      format      : 求解函数名
     *      initValue1  : 迭代初始值1
     *      initValue2  : 迭代初始值2, 1和2迭代初始值需不同但很近
     *      errMax      : 最大误差限
     *      MaxLoop     : 最大迭代次数
 * @ return 返回值: 
 *          求解结果
 */ 
REAL QuasiNewtonSolver(std::function<REAL(REAL)> format, REAL initValue1, REAL initValue2, REAL a = 0, REAL errMax = 1E-12, int MaxLoop = 20)
{   
    REAL xk = initValue1;
    REAL xkk = initValue2;
    REAL fk = format(xk) - a;
    REAL fkk = format(xkk) - a;
    REAL xkkk = xkk - (fkk * (xkk - xk))/(fkk - fk);
    REAL err = abs(xkkk - xkk);
    int loop = 1;
    while (loop < MaxLoop && err > errMax)
    {
        xk = xkk;
        xkk = xkkk;
        fk = fkk;
        fkk = format(xkkk) - a;
        xkkk = xkk - (fkk * (xkk - xk))/(fkk - fk);
        err = abs(xkkk - xkk);
        loop ++;
    }
    return xkkk;
}

#endif