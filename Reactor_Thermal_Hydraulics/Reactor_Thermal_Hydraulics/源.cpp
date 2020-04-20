#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define hl(p) (867.27 * p / 1000000 + 334.81) //饱和水焓hl的计算式
#define hll(p) (310.67 * p / 1000000 + 2645.6) //饱和蒸汽焓hll的计算式
#define wl(p) (-156.39 * p / 1000000 + 973.69) //饱和水密度w1的计算式
#define wll(p) (5.386 * p / 1000000 + 0.0536) //饱和蒸汽密度w11的计算式

//函数调用声明
double heat_exchanger(double M);
double two_phase_length(double M, double Ps, int n);
double down(double u);
double LAMDA(double Re);
double single_phase_h(double G, double Ps, double Pout);

//常数
const double g = 9.802, PI = 3.14159;
const int array_size = 20; //数组大小，可手动调
const double bar = 101325.0; //1 个大气压,Pa
const double Mu_l = 282.167e-6, Mu_s = 12.276e-6, Rho = 958.345;//1 个大气压下, 饱和水的汽、液相动力粘度Mu(μ) - Pa * s，水的密度ρ(Rho) - kg / m ^ 3
const double Ju = 0.131 + 0.163 * pow(0.3 / 0.3, 3.5);//；Ju—90°弯管的局部损失系数ξ = 0.291
const double Ju_tube = 1.5; //换热器进出、口局部损失系数ξ=1.5

//总管道几何尺寸
const double H0 = 4, H1 = 4, H2 = 6.5, H3 = 0.5, H4 = 5, H5 = 1, H6 = (H0 + H1 + H2 - H3 - H4 - H5);//H0 水箱液位；H4 换热器高度（换热管长度）；H6 竖直上升管长度, m
const double L1 = 3, L2 = 7, L3 = 4, D = 0.5;//横管长度,m // 管道管径0.3m
const double S = PI * D*D / 4; //主管道横截面积, m^2

//换热器几何参数
const int N_tube = 70; //N 为传热管根数
const double D_tube = 0.070, d_tube = 0.065;//D_tube 为传热管外径，d_tube传热管内径, 单位：m
const double S_tube = PI * d_tube*d_tube / 4; //S_tube 单根传热管流通横截面积, m ^ 2
const double A_tube = N_tube * PI*D_tube * H4; //A_tube 传热管外壁面面积, m ^ 2

//全局变量：从heat_exchanger 带回函数返回值
double Q;//换热器功率-MW
double delta_P_heat_exchanger;//换热器压降-Pa

//全局变量：统计整个回路压降
double Pg_total;
double Pa_total;
double Pf_total;
double Pj_total;

int main()
{
	int i = 0, n = 15;
	double M_OK[200] = {};
	double M, Ps, G;
	double P_in = 0, P_out = 0, T_in = 0, T_out = 0;
	double h_to_flash = 0, two_phase_h = 0;
	for (M = 1; M < 200; M++){
		printf("\n\n\n");
		printf("**************************************************流量M = %lf**********************************************\n", M);
		printf("*********************************************************************************************************\n\n");
		Pg_total = 0; Pa_total = 0; Pf_total = 0; Pj_total = 0;//归零，用于统计各项压降
		G = M / S;
		printf("------------------------------------------------------下降段---------------------------------------------- - \n");
		T_in = 100;
		P_in = down(G);
		printf("-----------------------------------------------------------------------------------------------------------\n\n");
		printf("------------------------------------------------------换热器---------------------------------------------- - \n");
		T_out = heat_exchanger(M); 
		if (T_out == -1) continue;
			P_out = P_in - delta_P_heat_exchanger;
		Ps = (T_out - 80.108)*1.0e6 / 204.81; //经验关系式：已知出口温度T_out，计算对应的饱和压力Ps
		if (P_out < Ps) {
			printf("ERROR: P_out < Ps 在换热器内就已发送闪蒸！！！\n"); 
			continue; 
		}
		printf("换热器出口温度：T_out=%.1lf ℃\n", T_out);
		printf("出口温度对应的饱和压力：Ps =%.0lf Pa\n", Ps);
		printf("换热器压降： Δ P =%.0lf Pa\n",delta_P_heat_exchanger);
		printf("换热器出口压力：P_out=%.0lf Pa\n", P_out);
		printf("换热器功率：Q =%.2lfMW\n", Q*1.0e-6);
		printf("-----------------------------------------------------------------------------------------------------------\n\n");
		printf("-----------------------------------------------------上升段（单相）---------------------------------------- - \n");
		h_to_flash = single_phase_h(G, Ps, P_out); if (h_to_flash == -1) continue;
		printf("换热器出口到闪蒸点的距离=%.2lf m\n", h_to_flash);
		printf("-----------------------------------------------------------------------------------------------------------\n\n");
		printf("-----------------------------------------------------上升段（两相）----------------------------------------\n");
		two_phase_h = two_phase_length(G, Ps, n);
		printf("两相段总长=%.2lfm\n", two_phase_h);
		printf("-----------------------------------------------------------------------------------------------------------\n\n");
		printf("-----------------------------------------------------------------------------------------------------------\n");
		printf("整个回路的重位压降：Pg_total =%-6.lf Pa\n", Pg_total);
		printf("整个回路的加速压降：Pa_total =%-6.lf Pa\n", Pa_total);
		printf("整个回路的摩阻压降：Pf_total =%-6.lf Pa\n", Pf_total);
		printf("整个回路的局部压降：Pj_total =%-6.lf Pa\n", Pj_total);
		printf("整个回路的总压降：P_total =%-6.lf Pa\n\n", Pg_total + Pa_total + Pf_total + Pj_total);
		printf("质量流量：M=%.2lf Kg/s\n", M);
		printf("流速：u=%.2lf m/s\n", M / Rho / S);
		printf("换热器功率：Q=%.2lfMW\n", Q*1.0e-6);
		printf("**********************************************************************************************************\n\n");
		if ((fabs((two_phase_h + h_to_flash) - (H5 + H6)) < 0.1) && (fabs(Pg_total + Pa_total + Pf_total + Pj_total) < 100)){
			i++; 
			M_OK[i] = M;
		}
		//break;
		//else if (((two_phase_h + h_to_flash) - (H5 + H6)) > 0.1)
		//if ((Pg_total + Pa_total + Pf_total + Pj_total > 50))
		//{ M = M + 1;}
		//else M = M - 1;
	}
	printf("******************************迭代结束******************************\n");
	if (i == 0) printf("没有合适的流量值！！！\n\n\n\n\n\n\n");
	for (; i >= 1; i--)
		printf("M=%.lf kg/s \t u=%.3lf m/s \n\n\n", M_OK[i],M_OK[i] / Rho / S);
	getchar();
	return 0;
}

//-------------------------------------------------------------------------------------------------
//函数功能：计算下降段（从水箱水面到换热器入口）总压降Pdown，得到换热器入口处压力值P_in，并返回该值
//注：
double down(double G)
{
	double Pdown, Pf, Pg, Pj; //Pdown-下降段出口压力值；Pf, Pg, Pj-摩阻压降、重位压降、局部压降
	double Re, lmd; //lmd-沿程损失系数λ
	Re = G * D / Mu_l;
	lmd = LAMDA(Re);
	Pg = Rho * g*(H0 + H1 + H2 - H3); //下降段重位压降
	Pf = lmd * (H1 + H2 + H3 + L1 + L2) / D * G*G / 2 / Rho; //下降段摩擦压降
	Pj = (4 * Ju + 0.5)*G*G / 2 / Rho; //下降段局部压降：4 个90°弯头ξ = 0.291 + 水箱进口1 个突缩ξ = 0.5
	Pdown = bar + Pg - Pf - Pj;
	//printf("下降段重位压降：Pg=%.0lf Pa\n", Pg);
	//printf("下降段摩擦压降：Pf=%.0lf Pa\n", Pf);
	//printf("下降段局部压降：Pj=%.0lf Pa\n", Pj);
	Pg_total = Pg_total + Pg;
	Pa_total = Pa_total + 0;
	Pf_total = Pf_total - Pf;
	Pj_total = Pj_total - Pj;
	printf("局部压降=%.lf \t 重位压降=%.lf \t 摩阻压降=%.lf \t 下降段出口处压力 = %.lf\n", -Pj, Pg, -Pf, (bar + Pg_total + Pa_total + Pf_total + Pj_total));
	return Pdown;
}

//函数功能：输入流量M，输出：换热器出口温度T_out、换热器换热功率Q、换热器压降P_hex
//
double heat_exchanger(double M)
{
	const double Ws = 0.54, p_mix = 430000.0, t_fo = 128.0;//Ws 为安全壳内蒸汽的质量份额，p 为混合气体压力，t_fo 为安全壳内混合气体温度，tw为传热管外壁面温度
	const double lamda_1 = 16.6, ts = 128.47, t_in = 100; //lamda_1 铬镍钢在100 摄氏度的导热系数；ts 安全壳内的水蒸气分压对应的饱和温度；t_in换热器入口温度
	const double Pr = 1.665, lamda_2 = 0.679, Cp = 4222.46, miu = 267.9e-6;
	/*按定性温度=105℃查得的：Pr 饱和水数、lamda_2 水的导热系数W/(m*K)、Cp 水的比热容(J/(kg*K))、miu 水的动力黏度为(Pa*s)*/
	double tw0 = 120, tw1; //传热管外壁面温度:：tw0 储存上次计算的值；tw1 储存新计算的值
	double t_out_0 = 108, t_out_1;//换热器出口温度：t_out_0 储存上次计算的值；t_out_1 储存新计算的值
	double ho, hi;//ho 传热管内侧传热系数；hi 管外侧传热系数
	double Re_tube, Nu;//管内侧传热系数：雷诺数、努塞尔数、普朗特数
	double R, K; //K 总传热系数;K = 1 / R
	double Q0, Q1; //Q0 = M*Cp*Δt(out-in)；Q1 = K*A*Δt_ln
	double delta_t;//对数平均温差
	do{
		t_out_1 = t_out_0;//t_out_0 储存上次计算的值，开始新一次计算ll
		//壁温迭代-------------------------------------------
		do{
			tw1 = tw0; //tw0 储存上次计算的值，开始新一次计算
			Q0 = M * Cp*(t_out_1 - t_in); //Q0=M*Cp*Δt(out-in)
			if (tw1 < t_fo) ho = 55.63*pow(Ws, 2.344)*pow(p_mix, 0.252)*pow(t_fo - tw1, 0.307);//计算管外侧换热系数ho
			else {
				printf("ERROR:壁温迭代出错,壁温>混合气体温度!!!\n");
				return -1;
			}
			tw0 = t_fo - Q0 / ho / A_tube;//tw1 储存新计算的值
		} while (fabs(tw1 - tw0) / tw1 > 0.001);
		//计算管外侧换热系数ho---------------------------------
		ho = 55.63*pow(Ws, 2.344)*pow(p_mix, 0.252)*pow(t_fo - tw1, 0.307);
		//计算管内侧换热系数hi----------------------------------
		Re_tube = 4 * M / (PI*d_tube*miu) / N_tube;
		Nu = 0.023*pow(Re_tube, 0.8)*pow(Pr, 0.4);
		hi = Nu * lamda_2 / d_tube;
		// 计算总传热系数---------------------------------------------------------------------------------------
		R = 1 / hi * D_tube / d_tube + 0.5*D_tube / lamda_1 * log(D_tube / d_tube) + 1 / ho; K = 1 / R;
		// 计算对数平均温差-------------------------------------------------------------------------------------
		if ((ts > t_out_1) && (t_out_1 > t_in)) 
			delta_t = (t_out_1 - t_in) / log((ts - t_in) / (ts - t_out_1));
		else {
				printf("ERROR:出口温度迭代出错!!! \t t_out_1=%.2lf℃ > 安全壳内温度128.47℃ 该流量不合适\n", t_out_1); 
				return -1; 
		}
		// 由K*A* Δ t_ln 计算新的传热量----------------------------------------
		Q1 = A_tube * K*delta_t;
		t_out_0 = Q1 / M / Cp + t_in;
	} while (fabs(t_out_1 - t_out_0) / t_out_0 > 0.001);
	//printf("壁温：tw=%.1lf℃\n", tw1);
	//printf("换热器出口温度：T_out=%.1lf℃\n", t_out_1);
	//printf("换热器功率：Q=%.2lfMW\n", Q1*1.0e-6);
	// 换热器压降计算-------------------------------------------------------------------------------------------
	double delta_P_g, delta_P_f, delta_P_j, lamda_tube; //P_hex 换热器总压降；lamda_tube 传热管内沿程损失系数
	double u_tube = M / (Rho*S_tube*N_tube); //换热管内的流速
	lamda_tube = LAMDA(Re_tube);
	delta_P_g = Rho * g*H4; //计算换热器重位压降
	delta_P_f = N_tube * lamda_tube*H4 / d_tube * Rho*u_tube*u_tube / 2;
	//计算摩擦压降
	delta_P_j = 2 * Ju_tube*Rho* (M / Rho / S)*(M / Rho / S) / 2;
	//进出口局部压降
	//printf("delta_P_g=%f\n", delta_P_g);
	//printf("delta_P_f=%f\n", delta_P_f);
	//printf("delta_P_j=%f\n", delta_P_j);
	//用全局变量带回多个返回值------
	delta_P_heat_exchanger = delta_P_g + delta_P_f + delta_P_j;
	Q = Q1;
	Pg_total = Pg_total - delta_P_g;
	Pa_total = Pa_total + 0;
	Pf_total = Pf_total - delta_P_f;
	Pj_total = Pj_total - delta_P_j;
	printf("局部压降=%.lf \t 重位压降=%.lf \t 摩阻压降=%.lf \t 换热器出口处压力 = %.lf\n", -delta_P_j, -delta_P_g, -delta_P_f, (bar + Pg_total + Pa_total + Pf_total + Pj_total));
	return t_out_1;
}

//----------------------------------------------
//函数功能：求换热器出口到闪蒸点的距离h
//注：
double single_phase_h(double G, double Ps, double Pout){
	double Pg, Pf, Pj, h = 0;//h 闪蒸点到竖直上升段起点的高度差
	double Re = G * D / Mu_l;
	double lmd = LAMDA(Re);
	double ppp;
	ppp = Pout - Rho * g*(H5)-2 * Ju*G*G / 2 / wl(Pout) - lmd * (H5 + L3) / D * G*G / 2 / wl(Pout);///ppp 是竖直上升段起点压力值，用来判断闪蒸能否发生在竖直上升段
	if (ppp < Ps) {
		printf("ERROR:进入竖直管段之前就已发生闪蒸！！！\n");
		return -1;
	}
	do{
		h += 0.01;
		Pg = Rho * g*(H5 + h); //换热器出口至闪蒸点段"重位压降"
		Pf = lmd * (H5 + L3 + h) / D * G*G / 2 / wl(Pout); //换热器出口至闪蒸点段"摩擦压降"
		Pj = 2 * Ju*G*G / 2 / wl(Pout); //换热器出口至闪蒸点段"局部压降"——2 个弯头
	} while (((Pout - Ps) - (Pg + Pf + Pj)) > 0);
	Pg_total = Pg_total - Pg;
	Pa_total = Pa_total + 0;
	Pf_total = Pf_total - Pf;
	Pj_total = Pj_total - Pj;
	printf("局部压降=%.lf \t 重位压降=%.lf \t 摩阻压降=%.lf \t 闪蒸点处压力 = %.lf\n", -Pj, -Pg, -Pf, (bar + Pg_total + Pa_total + Pf_total + Pj_total) );
	return (h + H5);
}

//------------------------------------------------------------------------------------------------
//函数功能：从闪蒸点开始，计算两相段理论长度L，并返回该值
//注：
double two_phase_length(double G, double Ps, int n){
	int i;
	double h0 = hl(Ps); //闪蒸点处的焓，即饱和压力Ps 下对应饱和液的焓值
	double dp = (Ps - bar) / n;//每个分段的压降为dp
	double p[array_size] = {}; p[0] = Ps; p[n] = bar;
	double x[array_size] = {}; x[0] = 0;
	double length[array_size] = {}, Total_length = 0;
	double Re_lo, lmd_lo, lo; //Re_lo 全液相的雷诺数；lmd_lo 全液相的摩阻系数；lo - 全液相折算系数Φ_lo
	double dpg, dpf, dpf_lo, Pa; //dpg-重位压降梯度；dpf-两相摩阻压降梯度；dpf_lo - 全液相摩擦压降梯度；Pa - 加速压降
	double pg_t = 0, pf_t = 0, pa_t = 0;//两相段的重位压降、摩阻压降、加速压降
	double Nu_m[array_size]; Nu_m[0] = 1 / wl(Ps); //Nu_m—流动比容ν_m, 是流动密度ρm 的倒数; 1 - 进口，2 - 出口；中间变量
	printf("每一段的压降dp=%6.1lf Pa\n", dp);
	for (i = 1; i <= n; i++){
		p[i] = p[0] - i * dp;
		x[i] = (h0 - hl(p[i])) / (hll(p[i]) - hl(p[i]));
		Nu_m[i] = x[i] / wll(p[i]) + (1 - x[i]) / wl(p[i]);
		//Nu_m_2 = x[i] / 0.590 + (1 - x[i]) / Rho;
		//Nu_m_1 = x[i - 1] / 0.590 + (1 - x[i - 1]) / Rho;
		//printf("第%d 段Nu_m_2=%.9lf\n", i, 1/Nu_m_2);
		//printf("第%d 段Nu_m_1=%.9lf\n", i,1/ Nu_m_1);
		//第i 段内的加速压降
		Pa = G * G*(Nu_m[i] - Nu_m[i - 1]);
		//第i 段内的重位压降梯度
		dpg = (1 / Nu_m[i] + 1 / Nu_m[i - 1]) / 2 * g; // (1 / Nu_m[i] + 1 / Nu_m[i - 1]) / 2 是进出口平均密度
		//dpg = Rho*g ;
		//第i 段内的两相摩阻压降梯度
		Re_lo = G * D / Mu_l; lmd_lo = LAMDA(Re_lo);
		dpf_lo = lmd_lo * G*G / D / 2 / Rho; //全液相摩擦压降梯度
		lo = ((wl(p[i]) / wll(p[i]) - 1)*x[i] + 1) * pow(((Mu_l / Mu_s	- 1) * x[i] + 1), -0.25);
		dpf = dpf_lo * lo;
		length[i] = (dp - Pa) / (dpg + dpf);
		Total_length = Total_length + length[i];
		//if (Total_length > two_phase_length_actual) break;
		pa_t += Pa;
		pg_t += dpg * length[i];
		pf_t += dpf * length[i];
		printf("第%2d 段：加速压降=%-6.lf \t 重位压降=%-6.lf \t 摩阻压降 = %-6.lf\n", i, -Pa, -dpg*length[i], -dpf*length[i]);
	}
	Pg_total = Pg_total - pg_t;
	Pa_total = Pa_total - pa_t;
	Pf_total = Pf_total - pf_t;
	Pj_total = Pj_total + 0;
	//printf("两相段总长=%.3lfm\n", Total_length);
	for (i = 0; i <= n; i++)
		printf("x[%2d]=%.10lf \t p[%2d]=%-6.f Pa \t Nu_m[%2d]=%.10lf \t Rho_m[% 2d] = %4.2lf \t length[% 2d] = %.2lf \n", i, x[i], i, p[i], i, Nu_m[i],	i, 1 / Nu_m[i], i, length[i]);
	printf("加速压降=%-6.lf \t 重位压降=%-6.lf \t 摩阻压降=%-6.lf \t 两相段出口处压力 = %-6.lf\n", -pa_t, -pg_t, -pf_t, (bar + Pg_total + Pa_total + Pf_total + Pj_total));
	return Total_length;
}

//-----------------------------------------
//函数功能：根据Re 数求沿程损失系数λ(lamda)
//注：
double LAMDA(double Re)
{
	double lamda;
	if (Re > 1 && Re < 2000) lamda = 64 / Re;
	else if (Re > 2000 && Re < 3.0E4) lamda = 0.3164*pow(Re, -0.25);
	else if (Re > 3.0E4) lamda = 0.184 * pow(Re, -0.2);
	else { 
		printf("雷诺数计算出错\n"); 
		exit(0); 
	}
	return lamda;
}