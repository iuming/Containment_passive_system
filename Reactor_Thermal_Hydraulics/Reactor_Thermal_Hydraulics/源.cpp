#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define hl(p) (867.27 * p / 1000000 + 334.81) //����ˮ��hl�ļ���ʽ
#define hll(p) (310.67 * p / 1000000 + 2645.6) //����������hll�ļ���ʽ
#define wl(p) (-156.39 * p / 1000000 + 973.69) //����ˮ�ܶ�w1�ļ���ʽ
#define wll(p) (5.386 * p / 1000000 + 0.0536) //���������ܶ�w11�ļ���ʽ

//������������
double heat_exchanger(double M);
double two_phase_length(double M, double Ps, int n);
double down(double u);
double LAMDA(double Re);
double single_phase_h(double G, double Ps, double Pout);

//����
const double g = 9.802, PI = 3.14159;
const int array_size = 20; //�����С�����ֶ���
const double bar = 101325.0; //1 ������ѹ,Pa
const double Mu_l = 282.167e-6, Mu_s = 12.276e-6, Rho = 958.345;//1 ������ѹ��, ����ˮ������Һ�ද��ճ��Mu(��) - Pa * s��ˮ���ܶȦ�(Rho) - kg / m ^ 3
const double Ju = 0.131 + 0.163 * pow(0.3 / 0.3, 3.5);//��Ju��90����ܵľֲ���ʧϵ���� = 0.291
const double Ju_tube = 1.5; //�������������ھֲ���ʧϵ����=1.5

//�ܹܵ����γߴ�
const double H0 = 4, H1 = 4, H2 = 6.5, H3 = 0.5, H4 = 5, H5 = 1, H6 = (H0 + H1 + H2 - H3 - H4 - H5);//H0 ˮ��Һλ��H4 �������߶ȣ����ȹܳ��ȣ���H6 ��ֱ�����ܳ���, m
const double L1 = 3, L2 = 7, L3 = 4, D = 0.5;//��ܳ���,m // �ܵ��ܾ�0.3m
const double S = PI * D*D / 4; //���ܵ�������, m^2

//���������β���
const int N_tube = 70; //N Ϊ���ȹܸ���
const double D_tube = 0.070, d_tube = 0.065;//D_tube Ϊ���ȹ��⾶��d_tube���ȹ��ھ�, ��λ��m
const double S_tube = PI * d_tube*d_tube / 4; //S_tube �������ȹ���ͨ������, m ^ 2
const double A_tube = N_tube * PI*D_tube * H4; //A_tube ���ȹ���������, m ^ 2

//ȫ�ֱ�������heat_exchanger ���غ�������ֵ
double Q;//����������-MW
double delta_P_heat_exchanger;//������ѹ��-Pa

//ȫ�ֱ�����ͳ��������·ѹ��
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
		printf("**************************************************����M = %lf**********************************************\n", M);
		printf("*********************************************************************************************************\n\n");
		Pg_total = 0; Pa_total = 0; Pf_total = 0; Pj_total = 0;//���㣬����ͳ�Ƹ���ѹ��
		G = M / S;
		printf("------------------------------------------------------�½���---------------------------------------------- - \n");
		T_in = 100;
		P_in = down(G);
		printf("-----------------------------------------------------------------------------------------------------------\n\n");
		printf("------------------------------------------------------������---------------------------------------------- - \n");
		T_out = heat_exchanger(M); 
		if (T_out == -1) continue;
			P_out = P_in - delta_P_heat_exchanger;
		Ps = (T_out - 80.108)*1.0e6 / 204.81; //�����ϵʽ����֪�����¶�T_out�������Ӧ�ı���ѹ��Ps
		if (P_out < Ps) {
			printf("ERROR: P_out < Ps �ڻ������ھ��ѷ�������������\n"); 
			continue; 
		}
		printf("�����������¶ȣ�T_out=%.1lf ��\n", T_out);
		printf("�����¶ȶ�Ӧ�ı���ѹ����Ps =%.0lf Pa\n", Ps);
		printf("������ѹ���� �� P =%.0lf Pa\n",delta_P_heat_exchanger);
		printf("����������ѹ����P_out=%.0lf Pa\n", P_out);
		printf("���������ʣ�Q =%.2lfMW\n", Q*1.0e-6);
		printf("-----------------------------------------------------------------------------------------------------------\n\n");
		printf("-----------------------------------------------------�����Σ����ࣩ---------------------------------------- - \n");
		h_to_flash = single_phase_h(G, Ps, P_out); if (h_to_flash == -1) continue;
		printf("���������ڵ�������ľ���=%.2lf m\n", h_to_flash);
		printf("-----------------------------------------------------------------------------------------------------------\n\n");
		printf("-----------------------------------------------------�����Σ����ࣩ----------------------------------------\n");
		two_phase_h = two_phase_length(G, Ps, n);
		printf("������ܳ�=%.2lfm\n", two_phase_h);
		printf("-----------------------------------------------------------------------------------------------------------\n\n");
		printf("-----------------------------------------------------------------------------------------------------------\n");
		printf("������·����λѹ����Pg_total =%-6.lf Pa\n", Pg_total);
		printf("������·�ļ���ѹ����Pa_total =%-6.lf Pa\n", Pa_total);
		printf("������·��Ħ��ѹ����Pf_total =%-6.lf Pa\n", Pf_total);
		printf("������·�ľֲ�ѹ����Pj_total =%-6.lf Pa\n", Pj_total);
		printf("������·����ѹ����P_total =%-6.lf Pa\n\n", Pg_total + Pa_total + Pf_total + Pj_total);
		printf("����������M=%.2lf Kg/s\n", M);
		printf("���٣�u=%.2lf m/s\n", M / Rho / S);
		printf("���������ʣ�Q=%.2lfMW\n", Q*1.0e-6);
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
	printf("******************************��������******************************\n");
	if (i == 0) printf("û�к��ʵ�����ֵ������\n\n\n\n\n\n\n");
	for (; i >= 1; i--)
		printf("M=%.lf kg/s \t u=%.3lf m/s \n\n\n", M_OK[i],M_OK[i] / Rho / S);
	getchar();
	return 0;
}

//-------------------------------------------------------------------------------------------------
//�������ܣ������½��Σ���ˮ��ˮ�浽��������ڣ���ѹ��Pdown���õ���������ڴ�ѹ��ֵP_in�������ظ�ֵ
//ע��
double down(double G)
{
	double Pdown, Pf, Pg, Pj; //Pdown-�½��γ���ѹ��ֵ��Pf, Pg, Pj-Ħ��ѹ������λѹ�����ֲ�ѹ��
	double Re, lmd; //lmd-�س���ʧϵ����
	Re = G * D / Mu_l;
	lmd = LAMDA(Re);
	Pg = Rho * g*(H0 + H1 + H2 - H3); //�½�����λѹ��
	Pf = lmd * (H1 + H2 + H3 + L1 + L2) / D * G*G / 2 / Rho; //�½���Ħ��ѹ��
	Pj = (4 * Ju + 0.5)*G*G / 2 / Rho; //�½��ξֲ�ѹ����4 ��90����ͷ�� = 0.291 + ˮ�����1 ��ͻ���� = 0.5
	Pdown = bar + Pg - Pf - Pj;
	//printf("�½�����λѹ����Pg=%.0lf Pa\n", Pg);
	//printf("�½���Ħ��ѹ����Pf=%.0lf Pa\n", Pf);
	//printf("�½��ξֲ�ѹ����Pj=%.0lf Pa\n", Pj);
	Pg_total = Pg_total + Pg;
	Pa_total = Pa_total + 0;
	Pf_total = Pf_total - Pf;
	Pj_total = Pj_total - Pj;
	printf("�ֲ�ѹ��=%.lf \t ��λѹ��=%.lf \t Ħ��ѹ��=%.lf \t �½��γ��ڴ�ѹ�� = %.lf\n", -Pj, Pg, -Pf, (bar + Pg_total + Pa_total + Pf_total + Pj_total));
	return Pdown;
}

//�������ܣ���������M������������������¶�T_out�����������ȹ���Q��������ѹ��P_hex
//
double heat_exchanger(double M)
{
	const double Ws = 0.54, p_mix = 430000.0, t_fo = 128.0;//Ws Ϊ��ȫ���������������ݶp Ϊ�������ѹ����t_fo Ϊ��ȫ���ڻ�������¶ȣ�twΪ���ȹ�������¶�
	const double lamda_1 = 16.6, ts = 128.47, t_in = 100; //lamda_1 ��������100 ���϶ȵĵ���ϵ����ts ��ȫ���ڵ�ˮ������ѹ��Ӧ�ı����¶ȣ�t_in����������¶�
	const double Pr = 1.665, lamda_2 = 0.679, Cp = 4222.46, miu = 267.9e-6;
	/*�������¶�=105���õģ�Pr ����ˮ����lamda_2 ˮ�ĵ���ϵ��W/(m*K)��Cp ˮ�ı�����(J/(kg*K))��miu ˮ�Ķ�����Ϊ(Pa*s)*/
	double tw0 = 120, tw1; //���ȹ�������¶�:��tw0 �����ϴμ����ֵ��tw1 �����¼����ֵ
	double t_out_0 = 108, t_out_1;//�����������¶ȣ�t_out_0 �����ϴμ����ֵ��t_out_1 �����¼����ֵ
	double ho, hi;//ho ���ȹ��ڲഫ��ϵ����hi ����ഫ��ϵ��
	double Re_tube, Nu;//���ڲഫ��ϵ������ŵ����Ŭ����������������
	double R, K; //K �ܴ���ϵ��;K = 1 / R
	double Q0, Q1; //Q0 = M*Cp*��t(out-in)��Q1 = K*A*��t_ln
	double delta_t;//����ƽ���²�
	do{
		t_out_1 = t_out_0;//t_out_0 �����ϴμ����ֵ����ʼ��һ�μ���ll
		//���µ���-------------------------------------------
		do{
			tw1 = tw0; //tw0 �����ϴμ����ֵ����ʼ��һ�μ���
			Q0 = M * Cp*(t_out_1 - t_in); //Q0=M*Cp*��t(out-in)
			if (tw1 < t_fo) ho = 55.63*pow(Ws, 2.344)*pow(p_mix, 0.252)*pow(t_fo - tw1, 0.307);//�������໻��ϵ��ho
			else {
				printf("ERROR:���µ�������,����>��������¶�!!!\n");
				return -1;
			}
			tw0 = t_fo - Q0 / ho / A_tube;//tw1 �����¼����ֵ
		} while (fabs(tw1 - tw0) / tw1 > 0.001);
		//�������໻��ϵ��ho---------------------------------
		ho = 55.63*pow(Ws, 2.344)*pow(p_mix, 0.252)*pow(t_fo - tw1, 0.307);
		//������ڲ໻��ϵ��hi----------------------------------
		Re_tube = 4 * M / (PI*d_tube*miu) / N_tube;
		Nu = 0.023*pow(Re_tube, 0.8)*pow(Pr, 0.4);
		hi = Nu * lamda_2 / d_tube;
		// �����ܴ���ϵ��---------------------------------------------------------------------------------------
		R = 1 / hi * D_tube / d_tube + 0.5*D_tube / lamda_1 * log(D_tube / d_tube) + 1 / ho; K = 1 / R;
		// �������ƽ���²�-------------------------------------------------------------------------------------
		if ((ts > t_out_1) && (t_out_1 > t_in)) 
			delta_t = (t_out_1 - t_in) / log((ts - t_in) / (ts - t_out_1));
		else {
				printf("ERROR:�����¶ȵ�������!!! \t t_out_1=%.2lf�� > ��ȫ�����¶�128.47�� ������������\n", t_out_1); 
				return -1; 
		}
		// ��K*A* �� t_ln �����µĴ�����----------------------------------------
		Q1 = A_tube * K*delta_t;
		t_out_0 = Q1 / M / Cp + t_in;
	} while (fabs(t_out_1 - t_out_0) / t_out_0 > 0.001);
	//printf("���£�tw=%.1lf��\n", tw1);
	//printf("�����������¶ȣ�T_out=%.1lf��\n", t_out_1);
	//printf("���������ʣ�Q=%.2lfMW\n", Q1*1.0e-6);
	// ������ѹ������-------------------------------------------------------------------------------------------
	double delta_P_g, delta_P_f, delta_P_j, lamda_tube; //P_hex ��������ѹ����lamda_tube ���ȹ����س���ʧϵ��
	double u_tube = M / (Rho*S_tube*N_tube); //���ȹ��ڵ�����
	lamda_tube = LAMDA(Re_tube);
	delta_P_g = Rho * g*H4; //���㻻������λѹ��
	delta_P_f = N_tube * lamda_tube*H4 / d_tube * Rho*u_tube*u_tube / 2;
	//����Ħ��ѹ��
	delta_P_j = 2 * Ju_tube*Rho* (M / Rho / S)*(M / Rho / S) / 2;
	//�����ھֲ�ѹ��
	//printf("delta_P_g=%f\n", delta_P_g);
	//printf("delta_P_f=%f\n", delta_P_f);
	//printf("delta_P_j=%f\n", delta_P_j);
	//��ȫ�ֱ������ض������ֵ------
	delta_P_heat_exchanger = delta_P_g + delta_P_f + delta_P_j;
	Q = Q1;
	Pg_total = Pg_total - delta_P_g;
	Pa_total = Pa_total + 0;
	Pf_total = Pf_total - delta_P_f;
	Pj_total = Pj_total - delta_P_j;
	printf("�ֲ�ѹ��=%.lf \t ��λѹ��=%.lf \t Ħ��ѹ��=%.lf \t ���������ڴ�ѹ�� = %.lf\n", -delta_P_j, -delta_P_g, -delta_P_f, (bar + Pg_total + Pa_total + Pf_total + Pj_total));
	return t_out_1;
}

//----------------------------------------------
//�������ܣ����������ڵ�������ľ���h
//ע��
double single_phase_h(double G, double Ps, double Pout){
	double Pg, Pf, Pj, h = 0;//h �����㵽��ֱ���������ĸ߶Ȳ�
	double Re = G * D / Mu_l;
	double lmd = LAMDA(Re);
	double ppp;
	ppp = Pout - Rho * g*(H5)-2 * Ju*G*G / 2 / wl(Pout) - lmd * (H5 + L3) / D * G*G / 2 / wl(Pout);///ppp ����ֱ���������ѹ��ֵ�������ж������ܷ�������ֱ������
	if (ppp < Ps) {
		printf("ERROR:������ֱ�ܶ�֮ǰ���ѷ�������������\n");
		return -1;
	}
	do{
		h += 0.01;
		Pg = Rho * g*(H5 + h); //�������������������"��λѹ��"
		Pf = lmd * (H5 + L3 + h) / D * G*G / 2 / wl(Pout); //�������������������"Ħ��ѹ��"
		Pj = 2 * Ju*G*G / 2 / wl(Pout); //�������������������"�ֲ�ѹ��"����2 ����ͷ
	} while (((Pout - Ps) - (Pg + Pf + Pj)) > 0);
	Pg_total = Pg_total - Pg;
	Pa_total = Pa_total + 0;
	Pf_total = Pf_total - Pf;
	Pj_total = Pj_total - Pj;
	printf("�ֲ�ѹ��=%.lf \t ��λѹ��=%.lf \t Ħ��ѹ��=%.lf \t �����㴦ѹ�� = %.lf\n", -Pj, -Pg, -Pf, (bar + Pg_total + Pa_total + Pf_total + Pj_total) );
	return (h + H5);
}

//------------------------------------------------------------------------------------------------
//�������ܣ��������㿪ʼ��������������۳���L�������ظ�ֵ
//ע��
double two_phase_length(double G, double Ps, int n){
	int i;
	double h0 = hl(Ps); //�����㴦���ʣ�������ѹ��Ps �¶�Ӧ����Һ����ֵ
	double dp = (Ps - bar) / n;//ÿ���ֶε�ѹ��Ϊdp
	double p[array_size] = {}; p[0] = Ps; p[n] = bar;
	double x[array_size] = {}; x[0] = 0;
	double length[array_size] = {}, Total_length = 0;
	double Re_lo, lmd_lo, lo; //Re_lo ȫҺ�����ŵ����lmd_lo ȫҺ���Ħ��ϵ����lo - ȫҺ������ϵ����_lo
	double dpg, dpf, dpf_lo, Pa; //dpg-��λѹ���ݶȣ�dpf-����Ħ��ѹ���ݶȣ�dpf_lo - ȫҺ��Ħ��ѹ���ݶȣ�Pa - ����ѹ��
	double pg_t = 0, pf_t = 0, pa_t = 0;//����ε���λѹ����Ħ��ѹ��������ѹ��
	double Nu_m[array_size]; Nu_m[0] = 1 / wl(Ps); //Nu_m���������ݦ�_m, �������ܶȦ�m �ĵ���; 1 - ���ڣ�2 - ���ڣ��м����
	printf("ÿһ�ε�ѹ��dp=%6.1lf Pa\n", dp);
	for (i = 1; i <= n; i++){
		p[i] = p[0] - i * dp;
		x[i] = (h0 - hl(p[i])) / (hll(p[i]) - hl(p[i]));
		Nu_m[i] = x[i] / wll(p[i]) + (1 - x[i]) / wl(p[i]);
		//Nu_m_2 = x[i] / 0.590 + (1 - x[i]) / Rho;
		//Nu_m_1 = x[i - 1] / 0.590 + (1 - x[i - 1]) / Rho;
		//printf("��%d ��Nu_m_2=%.9lf\n", i, 1/Nu_m_2);
		//printf("��%d ��Nu_m_1=%.9lf\n", i,1/ Nu_m_1);
		//��i ���ڵļ���ѹ��
		Pa = G * G*(Nu_m[i] - Nu_m[i - 1]);
		//��i ���ڵ���λѹ���ݶ�
		dpg = (1 / Nu_m[i] + 1 / Nu_m[i - 1]) / 2 * g; // (1 / Nu_m[i] + 1 / Nu_m[i - 1]) / 2 �ǽ�����ƽ���ܶ�
		//dpg = Rho*g ;
		//��i ���ڵ�����Ħ��ѹ���ݶ�
		Re_lo = G * D / Mu_l; lmd_lo = LAMDA(Re_lo);
		dpf_lo = lmd_lo * G*G / D / 2 / Rho; //ȫҺ��Ħ��ѹ���ݶ�
		lo = ((wl(p[i]) / wll(p[i]) - 1)*x[i] + 1) * pow(((Mu_l / Mu_s	- 1) * x[i] + 1), -0.25);
		dpf = dpf_lo * lo;
		length[i] = (dp - Pa) / (dpg + dpf);
		Total_length = Total_length + length[i];
		//if (Total_length > two_phase_length_actual) break;
		pa_t += Pa;
		pg_t += dpg * length[i];
		pf_t += dpf * length[i];
		printf("��%2d �Σ�����ѹ��=%-6.lf \t ��λѹ��=%-6.lf \t Ħ��ѹ�� = %-6.lf\n", i, -Pa, -dpg*length[i], -dpf*length[i]);
	}
	Pg_total = Pg_total - pg_t;
	Pa_total = Pa_total - pa_t;
	Pf_total = Pf_total - pf_t;
	Pj_total = Pj_total + 0;
	//printf("������ܳ�=%.3lfm\n", Total_length);
	for (i = 0; i <= n; i++)
		printf("x[%2d]=%.10lf \t p[%2d]=%-6.f Pa \t Nu_m[%2d]=%.10lf \t Rho_m[% 2d] = %4.2lf \t length[% 2d] = %.2lf \n", i, x[i], i, p[i], i, Nu_m[i],	i, 1 / Nu_m[i], i, length[i]);
	printf("����ѹ��=%-6.lf \t ��λѹ��=%-6.lf \t Ħ��ѹ��=%-6.lf \t ����γ��ڴ�ѹ�� = %-6.lf\n", -pa_t, -pg_t, -pf_t, (bar + Pg_total + Pa_total + Pf_total + Pj_total));
	return Total_length;
}

//-----------------------------------------
//�������ܣ�����Re �����س���ʧϵ����(lamda)
//ע��
double LAMDA(double Re)
{
	double lamda;
	if (Re > 1 && Re < 2000) lamda = 64 / Re;
	else if (Re > 2000 && Re < 3.0E4) lamda = 0.3164*pow(Re, -0.25);
	else if (Re > 3.0E4) lamda = 0.184 * pow(Re, -0.2);
	else { 
		printf("��ŵ���������\n"); 
		exit(0); 
	}
	return lamda;
}