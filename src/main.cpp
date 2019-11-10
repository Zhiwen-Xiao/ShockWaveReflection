/*
 * @version: v0.1
 * @Descripttion: 
 * @Author: XiaoZhiwen
 * @Date: 2019-11-09 17:13:12
 */
#include "inc.hpp"

typedef long double REAL;

#define Pi 3.141592653589793115997963468544185161590576171875

int main()
{
     using namespace std;

     REAL Ma, P0, Rho0, T0, theta, gamma, beta1, beta2, beta3, thetaS, ma1, ma2, ma3, p1, p2, p3, t1, t2, t3, rho1, rho2, rho3;
     
     //读取文件，输出
     cout.setf(ios::fixed);
     ifstream fin("input.txt");
     ofstream fout("output.txt");
     fout.setf(ios::fixed);
     fout << "输入参数：" <<endl;
     const int LIMIT = 255;
     char input[LIMIT];
     fin.getline(input,256);
     fin.getline(input,256);
     Ma = atof(input);
     cout << "Ma \t = \t" << setprecision(15) << Ma << endl;
     fout << "来流(0区)马赫数：\t\t\t\t" << setprecision(6) << Ma << endl;
     fin.getline(input,256);
     fin.getline(input,256);
     P0 = atof(input);
     cout << "P0 \t = \t" << setprecision(15) << P0 << endl;
     fout << "来流(0区)压力(Pa)：\t\t\t\t" << setprecision(6) << P0 << endl;
     fin.getline(input,256);
     fin.getline(input,256);
     Rho0 = atof(input);
     cout << "Rho0 \t = \t" << setprecision(15) << Rho0 << endl;
     fout << "来流(0区)密度(kg/cum)：\t\t\t" << setprecision(6) << Rho0 << endl;
     fin.getline(input,256);
     fin.getline(input,256);
     T0 = atof(input);
     cout << "T0 \t = \t" << setprecision(15) << T0 << endl;
     fout << "来流(0区)温度(K)：\t\t\t\t" << setprecision(6) << T0 << endl;
     fin.getline(input,256);
     fin.getline(input,256);
     theta = atof(input);
     cout << "theta \t = \t" << setprecision(15) << theta << endl;
     fout << "楔角(deg)：\t\t\t\t\t\t" << setprecision(6) << theta << endl;
     theta = theta / 180 * Pi;
     fin.getline(input,256);
     fin.getline(input,256);
     gamma = atof(input);
     cout << "gamma \t = \t" << setprecision(15) << gamma << endl;
     fout << "气体比热比：\t\t\t\t\t" << setprecision(6) << gamma << endl;
     fout << endl;


     REAL machConditionTheta = MachConditionTheta(Ma, gamma);
     REAL vonNeumannConditionTheta = VonNeumannConditionTheta(Ma, gamma);
     if(vonNeumannConditionTheta > machConditionTheta) 
     {
          vonNeumannConditionTheta = machConditionTheta;
     }

     if (theta <= vonNeumannConditionTheta)
     {
          fout << "楔角小于给定条件下冯诺依曼给出的马赫反射压力平衡条件（ theta_w^N = " << vonNeumannConditionTheta / Pi * 180 << " )，只发生正规反射" << endl;

          beta1 = Beta(theta, Ma, 0, gamma);
          ma1 = Ma2(Ma, beta1, gamma);
          beta2 = Beta(theta, ma1, 0, gamma);
          p1 = P0 * P_Ratio(Ma, beta1, gamma);
          rho1 = Rho0 * Rho_Ratio(Ma, beta1, gamma);
          t1 = T0 * T_Ratio(Ma, beta1, gamma);
          p2 = p1 * P_Ratio(ma1, beta2, gamma);
          rho2 = rho1 * Rho_Ratio(ma1, beta2, gamma);
          t2 = t1 * T_Ratio(ma1, beta2, gamma);
          ma2 = Ma2(ma1, beta2, gamma);

          fout << "入射激波激波角(deg)：\t\t\t" << setprecision(6) << beta1 / Pi * 180 << endl;
          fout << "反射激波激波角(deg)：\t\t\t" << setprecision(6) << beta2 / Pi * 180 << endl;
          fout << "\n入射激波后(1区)参数" << endl;
          fout << "入射激波后(1区)马赫数：\t\t\t" << setprecision(6) << ma1 << endl;
          fout << "入射激波后(1区)压力(Pa)：\t\t" << setprecision(6) << p1 << endl;
          fout << "入射激波后(1区)密度(kg/cum)：\t" << setprecision(6) << rho1 << endl;
          fout << "入射激波后(1区)温度(K)：\t\t" << setprecision(6) << t1 << endl;
          fout << "\n反射激波后(2区)参数" << endl;
          fout << "反射激波后(2区)马赫数：\t\t\t" << setprecision(6) << ma2 << endl;
          fout << "反射激波后(2区)压力(Pa)：\t\t" << setprecision(6) << p2 << endl;
          fout << "反射激波后(2区)密度(kg/cum)：\t" << setprecision(6) << rho2 << endl;
          fout << "反射激波后(2区)温度(K)：\t\t" << setprecision(6) << t2 << endl;
     } else if (theta >= machConditionTheta)
     {
          fout << "楔角大于给定条件下马赫反射的脱体条件（ theta_w^D = " << machConditionTheta / Pi * 180 << " )，只发生马赫反射" << endl;

          beta1 = Beta(theta, Ma, 0, gamma);
          ma1 = Ma2(Ma, beta1, gamma);
          thetaS = ThetaS_M(theta, Ma, gamma);
          beta2 = Beta(theta - thetaS, ma1, 0, gamma);
          beta3 = Beta(thetaS, Ma, 1, gamma);
          p1 = P0 * P_Ratio(Ma, beta1, gamma);
          rho1 = Rho0 * Rho_Ratio(Ma, beta1, gamma);
          t1 = T0 * T_Ratio(Ma, beta1, gamma);
          p2 = p1 * P_Ratio(ma1, beta2, gamma);
          rho2 = rho1 * Rho_Ratio(ma1, beta2, gamma);
          t2 = t1 * T_Ratio(ma1, beta2, gamma);
          ma2 = Ma2(ma1, beta2, gamma);
          p3 = P0 * P_Ratio(Ma, beta3, gamma);
          rho3 = Rho0 * Rho_Ratio(Ma, beta3, gamma);
          t3 = T0 * T_Ratio(Ma, beta3, gamma);
          ma3 = Ma2(Ma, beta3, gamma);

          fout << "入射激波激波角(deg)：\t\t\t" << setprecision(6) << beta1 / Pi * 180 << endl;
          fout << "反射激波激波角(deg)：\t\t\t" << setprecision(6) << beta2 / Pi * 180 << endl;
          fout << "马赫杆激波角(deg)：\t\t\t\t" << setprecision(6) << beta3 / Pi * 180 << endl;
          fout << "滑移线下偏角(deg)：\t\t\t\t" << setprecision(6) << thetaS / Pi * 180 << endl;
          fout << "\n入射激波后(1区)参数" << endl;
          fout << "入射激波后(1区)马赫数：\t\t\t" << setprecision(6) << ma1 << endl;
          fout << "入射激波后(1区)压力(Pa)：\t\t" << setprecision(6) << p1 << endl;
          fout << "入射激波后(1区)密度(kg/cum)：\t" << setprecision(6) << rho1 << endl;
          fout << "入射激波后(1区)温度(K)：\t\t" << setprecision(6) << t1 << endl;
          fout << "\n反射激波后(2区)参数" << endl;
          fout << "反射激波后(2区)马赫数：\t\t\t" << setprecision(6) << ma2 << endl;
          fout << "反射激波后(2区)压力(Pa)：\t\t" << setprecision(6) << p2 << endl;
          fout << "反射激波后(2区)密度(kg/cum)：\t" << setprecision(6) << rho2 << endl;
          fout << "反射激波后(2区)温度(K)：\t\t" << setprecision(6) << t2 << endl;
          fout << "\n马赫杆后(3区)参数" << endl;
          fout << "马赫杆后(3区)马赫数： \t\t\t" << setprecision(6) << ma3 << endl;
          fout << "马赫杆后(3区)压力(Pa)：\t\t\t" << setprecision(6) << p3 << endl;
          fout << "马赫杆后(3区)密度(kg/cum)：\t\t" << setprecision(6) << rho3 << endl;
          fout << "马赫杆后(3区)温度(K)：\t\t\t" << setprecision(6) << t3 << endl;
     }else
     {
          fout << "楔角大于给定条件下冯诺依曼给出的马赫反射压力平衡条件（ theta_w^N = " << vonNeumannConditionTheta / Pi * 180 << " ),而小于马赫反射的脱体条件（ theta_w^D = " << machConditionTheta / Pi * 180 << " )，故处于双解区" << endl;

          fout << "\n当楔角由小变大达到该条件，发生正规反射时：\n" << endl;

          beta1 = Beta(theta, Ma, 0, gamma);
          ma1 = Ma2(Ma, beta1, gamma);
          beta2 = Beta(theta, ma1, 0, gamma);
          p1 = P0 * P_Ratio(Ma, beta1, gamma);
          rho1 = Rho0 * Rho_Ratio(Ma, beta1, gamma);
          t1 = T0 * T_Ratio(Ma, beta1, gamma);
          p2 = p1 * P_Ratio(ma1, beta2, gamma);
          rho2 = rho1 * Rho_Ratio(ma1, beta2, gamma);
          t2 = t1 * T_Ratio(ma1, beta2, gamma);
          ma2 = Ma2(ma1, beta2, gamma);

          fout << "入射激波激波角(deg)：\t\t\t" << setprecision(6) << beta1 / Pi * 180 << endl;
          fout << "反射激波激波角(deg)：\t\t\t" << setprecision(6) << beta2 / Pi * 180 << endl;
          fout << "\n入射激波后(1区)参数" << endl;
          fout << "入射激波后(1区)马赫数：\t\t\t" << setprecision(6) << ma1 << endl;
          fout << "入射激波后(1区)压力(Pa)：\t\t" << setprecision(6) << p1 << endl;
          fout << "入射激波后(1区)密度(kg/cum)：\t" << setprecision(6) << rho1 << endl;
          fout << "入射激波后(1区)温度(K)：\t\t" << setprecision(6) << t1 << endl;
          fout << "\n反射激波后(2区)参数" << endl;
          fout << "反射激波后(2区)马赫数：\t\t\t" << setprecision(6) << ma2 << endl;
          fout << "反射激波后(2区)压力(Pa)：\t\t" << setprecision(6) << p2 << endl;
          fout << "反射激波后(2区)密度(kg/cum)：\t" << setprecision(6) << rho2 << endl;
          fout << "反射激波后(2区)温度(K)：\t\t" << setprecision(6) << t2 << endl;

          fout << "\n\n当楔角由大变小达到该条件，发生马赫反射时：\n" << endl;

          beta1 = Beta(theta, Ma, 0, gamma);
          ma1 = Ma2(Ma, beta1, gamma);
          thetaS = ThetaS_M(theta, Ma, gamma);
          beta2 = Beta(theta - thetaS, ma1, 0, gamma);
          beta3 = Beta(thetaS, Ma, 1, gamma);
          p1 = P0 * P_Ratio(Ma, beta1, gamma);
          rho1 = Rho0 * Rho_Ratio(Ma, beta1, gamma);
          t1 = T0 * T_Ratio(Ma, beta1, gamma);
          p2 = p1 * P_Ratio(ma1, beta2, gamma);
          rho2 = rho1 * Rho_Ratio(ma1, beta2, gamma);
          t2 = t1 * T_Ratio(ma1, beta2, gamma);
          ma2 = Ma2(ma1, beta2, gamma);
          p3 = P0 * P_Ratio(Ma, beta3, gamma);
          rho3 = Rho0 * Rho_Ratio(Ma, beta3, gamma);
          t3 = T0 * T_Ratio(Ma, beta3, gamma);
          ma3 = Ma2(Ma, beta3, gamma);

          fout << "入射激波激波角(deg)：\t\t\t" << setprecision(6) << beta1 / Pi * 180 << endl;
          fout << "反射激波激波角(deg)：\t\t\t" << setprecision(6) << beta2 / Pi * 180 << endl;
          fout << "马赫杆激波角(deg)：\t\t\t\t" << setprecision(6) << beta3 / Pi * 180 << endl;
          fout << "滑移线下偏角(deg)：\t\t\t\t" << setprecision(6) << thetaS / Pi * 180 << endl;
          fout << "\n入射激波后(1区)参数" << endl;
          fout << "入射激波后(1区)马赫数：\t\t\t" << setprecision(6) << ma1 << endl;
          fout << "入射激波后(1区)压力(Pa)：\t\t" << setprecision(6) << p1 << endl;
          fout << "入射激波后(1区)密度(kg/cum)：\t" << setprecision(6) << rho1 << endl;
          fout << "入射激波后(1区)温度(K)：\t\t" << setprecision(6) << t1 << endl;
          fout << "\n反射激波后(2区)参数" << endl;
          fout << "反射激波后(2区)马赫数：\t\t\t" << setprecision(6) << ma2 << endl;
          fout << "反射激波后(2区)压力(Pa)：\t\t" << setprecision(6) << p2 << endl;
          fout << "反射激波后(2区)密度(kg/cum)：\t" << setprecision(6) << rho2 << endl;
          fout << "反射激波后(2区)温度(K)：\t\t" << setprecision(6) << t2 << endl;
          fout << "\n马赫杆后(3区)参数" << endl;
          fout << "马赫杆后(3区)马赫数： \t\t\t" << setprecision(6) << ma3 << endl;
          fout << "马赫杆后(3区)压力(Pa)：\t\t\t" << setprecision(6) << p3 << endl;
          fout << "马赫杆后(3区)密度(kg/cum)：\t\t" << setprecision(6) << rho3 << endl;
          fout << "马赫杆后(3区)温度(K)：\t\t\t" << setprecision(6) << t3 << endl;
     }
     fin.close();
     fout.close();
     return 0;
}