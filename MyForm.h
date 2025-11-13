#pragma once
#include <math.h>
#include <tuple>
#include <iostream>
#include <vector>
const int p = 4; // р - порядок метода Рунге Кутта

/*	
*	Функция Runge_Kytta_4 - это функция для вычисления следующей точки численной траектории методом Рунге Кутта 4 порядка
*	Возвращает следующую точку численной траектории { x_n, v_n } 
*	double(*f)(double, double) - функция правой части дифференциального уравнения 
*	double h_n - шаг для изменения х
*	double x_n - значение х текущей точки 
*	double v_n - значение v текущей точки 
*/
std::pair<double, double> Runge_Kytta_4(double(*f)(double, double), double h_n, double x_n, double v_n) {

	double k1 = f(x_n, v_n);
	double k2 = f(x_n + h_n / 2.0, v_n + h_n / 2.0 * k1);
	double k3 = f(x_n + h_n / 2.0, v_n + h_n / 2.0 * k2);
	double k4 = f(x_n + h_n, v_n + h_n * k3);

	x_n = x_n + h_n;
	v_n = v_n + h_n * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
	return { x_n, v_n };
}
/*	
*	Функция Runge_Kytta_4_system - это функция для вычисления следующей точки численной траектории методом Рунге Кутта 4 порядка
*	для системы дифференциальных уравнений
*	Возвращает tuple, который хранит численное решение для переменно v1 и v2, {x_n, v1_n, v2_n}
*	std::pair<double, double>(*f)(double, double, double, double, double) - функции правой части системы дифференциальных уравнений
*	double h_n - шаг для изменения х
*	double x_n - значение х текущей точки 
*	double v1_n - значение v1 в текущей точке
*	double v2_n - значение v2 в текущей точке
*	double a, double b - коэффициенты системы
*/
std::tuple<double, double, double>Runge_Kytta_4_system(std::pair<double, double>(*f)(double, double, double, double, double), double h_n, double x_n, double v1_n, double v2_n, double a, double b) {
	double k1_u1 = f(x_n, v1_n, v2_n, a, b).first;
	double k1_u2 = f(x_n, v1_n, v2_n, a, b).second;

	double k2_u1 = f(x_n + h_n / 2.0, v1_n + h_n / 2.0 * k1_u1, v2_n + h_n / 2.0 * k1_u2, a, b).first;
	double k2_u2 = f(x_n + h_n / 2.0, v1_n + h_n / 2.0 * k1_u1, v2_n + h_n / 2.0 * k1_u2, a, b).second;

	double k3_u1 = f(x_n + h_n / 2.0, v1_n + h_n / 2.0 * k2_u1, v2_n + h_n / 2.0 * k2_u2, a, b).first;
	double k3_u2 = f(x_n + h_n / 2.0, v1_n + h_n / 2.0 * k2_u1, v2_n + h_n / 2.0 * k2_u2, a, b).second;

	double k4_u1 = f(x_n + h_n, v1_n + h_n * k3_u1, v2_n + h_n * k3_u2, a, b).first;
	double k4_u2 = f(x_n + h_n, v1_n + h_n * k3_u1, v2_n + h_n * k3_u2, a, b).second;

	x_n = x_n + h_n;
	v1_n = v1_n + h_n * (k1_u1 + 2.0 * k2_u1 + 2.0 * k3_u1 + k4_u1) / 6.0;
	v2_n = v2_n + h_n * (k1_u2 + 2.0 * k2_u2 + 2.0 * k3_u2 + k4_u2) / 6.0;

	return {x_n, v1_n, v2_n};
}
/*
*	Функция true_trajectory - истинное решение тестовой задачи 
*	Возвращает значение функции в точке х
*	double x - значение х текущей точки
*	double u0 - начальное условие 
*/
double true_trajectory(double x, double u0) {
	return u0 * exp(-2.5 * x);
}
/*
*	Функция test_function - функция, для которой строиться численная траектория, тестовая задача
*	Возвращает значение функции в точке
*/
double test_function(double x, double v) {
	return -2.5 * v;
}
/*
*	Функция function_1 - функция, для которой строиться численная траектория, задача 1
*	Возвращает значение функции в точке
*/
double function_1(double x, double v) {
	return (std::log(x + 1) / (pow(x, 2) + 1)) * pow(v, 2) + v - pow(v, 3) * sin(10 * x);
}
/*
*	Функция function_2 - система функция, для которой строиться численная траектория, задача 2
*	Возвращает значение функции в точке
*/
std::pair<double , double> function_2(double x, double u1, double u2, double a, double b) {
	double du1 = u2;
	double du2 = -a * pow(u2, 2) - b * sin(u1);

	return { du1, du2 };
}
/*
*	Функция S - функция нахождения контрольной величины
*	Возвращает значение контрольной величины 
*	double v_n - значение численной траектории, найденной с полным шагом
*	double v - значение численной траектории, найденной с помощью счета с половинным шагом
*/
double S(double v_n, double v) {
	return abs((v_n - v)) / (pow(2, p) - 1);
}
/*
*	Функция RK_4_OLP - функция, выпоняющая метод регулировки шага, для тестовой задачи и задачи 1
*	double(*f)(double, double) - функция правой части дифференциального уравнения
*	double x0 - значение х текущей точки 
*	double u0 - значение v текущей точки 
*	double h - шаг для изменения х
*	double e - точность локальной погрешности
*/
std::vector<double> RK_4_OLP(double(*f)(double, double), double x0, double u0, double h, double e) {
	std::pair<double, double> new_point_h;
	std::pair<double, double> new_point_2h;
	bool flag = true;
	double S_new = 0;
	double swich = 0;

	while (flag) {
		new_point_h = Runge_Kytta_4(f, h, x0, u0);
		new_point_2h = Runge_Kytta_4(f, h / 2.0, x0 + h/2.0, Runge_Kytta_4(f, h / 2.0, x0, u0).second);
		S_new = S(new_point_h.second, new_point_2h.second);

		if (S_new > e) {
			h /= 2.0;
			swich -= 1;
		} else
		if (S_new < e / pow(2, p + 1)) {
			flag = false;
			h *= 2.0;
			swich += 1;
		} else 
		if(S_new >= e / pow(2, p + 1) && S_new < e)
			flag = false;
	}
	std::vector<double> result = { new_point_h.first, new_point_h.second, new_point_2h.second, h, swich};
	return result;
}
/*
*	Функция RK_4_OLP_for_system - функция, выпоняющая метод регулировки шага для систем задачи 2
*	std::pair<double, double>(*f)(double, double, double, double, double) - функция правой части дифференциального уравнения
*	double x0 - значение х текущей точки 
*	double u0_1 - значение v1 в текущей точке
*	double u0_2 - значение v2 в текущей точке
*	double h - шаг для изменения х
*	double e - точность локальной погрешности
*	double a, double b - коэффициенты системы
*/
std::vector<double> RK_4_OLP_for_system(std::pair<double, double>(*f)(double, double, double, double, double), double x0, double u0_1, double u0_2, double h, double e, double a, double b) {
	bool flag = true;
	std::tuple<double, double, double> new_point_h;
	std::tuple<double, double, double> new_point_2h;
	double S1_new = 0;
	double S2_new = 0;
	double swich = 0;

	while (flag) {
		new_point_h = Runge_Kytta_4_system(f, h, x0, u0_1, u0_2, a, b);
		new_point_2h = Runge_Kytta_4_system(f, h / 2.0, x0 + h / 2.0, std::get<1>(Runge_Kytta_4_system(f, h / 2.0, x0, u0_1, u0_2, a, b)), std::get<2>(Runge_Kytta_4_system(f, h / 2.0, x0, u0_1, u0_2, a, b)), a, b);
		S1_new = S(std::get<1>(new_point_h), std::get<1>(new_point_2h));
		S2_new = S(std::get<2>(new_point_h), std::get<2>(new_point_2h));

		if (S1_new > e || S2_new > e) {
			h /= 2.0;
			swich -= 1;
		}
		else if (S1_new < e / pow(2, p + 1) && S2_new < e / pow(2, p + 1)) {
			flag = false;
			h *= 2.0;
			swich += 1;
		}
		else if ((S1_new >= e / pow(2, p + 1) && S1_new < e) || (S2_new >= e / pow(2, p + 1) && S2_new < e))
			flag = false;
	}
	std::vector<double> result = { std::get<0>(new_point_h), std::get<1>(new_point_h), std::get<2>(new_point_h), std::get<1>(new_point_2h), std::get<2>(new_point_2h), h, swich };
	return result;
}
namespace Graph {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;
	using namespace ZedGraph;

	/// <summary>
	/// Summary for MyForm
	/// </summary>
	public ref class MyForm : public System::Windows::Forms::Form
	{
	public:
		MyForm(void)
		{
			InitializeComponent();
			//
			//TODO: Add the constructor code here
			//
		}

	protected:
		/// <summary>
		/// Clean up any resources being used.
		/// </summary>
		~MyForm()
		{
			if (components)
			{
				delete components;
			}
		}
	private: ZedGraph::ZedGraphControl^ zedGraphControl1;
		   //private: ZedGraph::ZedGraphControl^  zedGraphControl2;
	private: System::Windows::Forms::Button^ button1;
	private: System::Windows::Forms::Button^ button3;
	private: System::Windows::Forms::DataGridView^ dataGridView1;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ i;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ x_i;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ v_i;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ v2_i;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ v_i_v2_i;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ OLP;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ h_i;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ C1;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ C2;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ u_i;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ u_i_v_i;
	private: System::Windows::Forms::Label^ label1;
	private: System::Windows::Forms::TextBox^ textBox1;
	private: System::Windows::Forms::Label^ label2;
	private: System::Windows::Forms::TextBox^ textBox2;
	private: System::Windows::Forms::Label^ label3;
	private: System::Windows::Forms::TextBox^ textBox3;
	private: System::Windows::Forms::Button^ button2;
	private: System::Windows::Forms::TextBox^ textBox6;
	private: System::Windows::Forms::Label^ label6;
	private: System::Windows::Forms::TextBox^ textBox7;
	private: System::Windows::Forms::Label^ label7;
	private: System::Windows::Forms::TextBox^ textBox4;
	private: System::Windows::Forms::Label^ label4;
	private: System::Windows::Forms::Button^ button4;
	private: System::Windows::Forms::Label^ label5;
	private: System::Windows::Forms::Label^ label8;
	private: System::Windows::Forms::Label^ label9;
	private: System::Windows::Forms::Button^ button5;
	private: System::Windows::Forms::Button^ button6;
	private: System::Windows::Forms::Button^ button7;
	private: System::Windows::Forms::Label^ label10;
	private: System::Windows::Forms::Label^ label11;
	private: System::Windows::Forms::TextBox^ textBox5;
	private: System::Windows::Forms::TextBox^ textBox8;
	private: System::Windows::Forms::Label^ label12;
	private: System::Windows::Forms::Label^ label13;
	private: System::Windows::Forms::TextBox^ textBox9;
	private: System::Windows::Forms::TextBox^ textBox10;
	private: System::Windows::Forms::TextBox^ textBox11;
	private: System::Windows::Forms::Label^ label14;
	private: System::Windows::Forms::TextBox^ textBox12;
	private: System::Windows::Forms::TextBox^ textBox13;
	private: System::Windows::Forms::TextBox^ textBox14;
	private: System::Windows::Forms::Label^ label15;
	private: System::Windows::Forms::Label^ label16;
	private: System::Windows::Forms::Label^ label17;
	private: System::Windows::Forms::Label^ label18;
	private: System::Windows::Forms::Label^ label19;
	private: System::Windows::Forms::TextBox^ textBox15;
	private: System::Windows::Forms::TextBox^ textBox16;
	private: System::Windows::Forms::TextBox^ textBox17;
	private: System::Windows::Forms::TextBox^ textBox18;
	private: System::Windows::Forms::Label^ label20;
	private: System::Windows::Forms::Label^ label21;
	private: System::Windows::Forms::Label^ label22;
	private: System::Windows::Forms::Label^ label23;
	private: System::Windows::Forms::TextBox^ textBox19;
	private: System::Windows::Forms::Label^ label24;
	private: System::Windows::Forms::Label^ label25;
	private: System::Windows::Forms::Label^ label26;
	private: System::Windows::Forms::TextBox^ textBox20;
	private: System::Windows::Forms::TextBox^ textBox21;
	private: System::Windows::Forms::Label^ label27;
	private: System::Windows::Forms::TextBox^ textBox22;
	private: System::Windows::Forms::Label^ label28;
	private: System::Windows::Forms::TextBox^ textBox23;
	protected:
	private: System::ComponentModel::IContainer^ components;

	private:
		/// <summary>
		/// Required designer variable.
		/// </summary>


#pragma region Windows Form Designer generated code
		/// <summary>
		/// Required method for Designer support - do not modify
		/// the contents of this method with the code editor.
		/// </summary>
		void InitializeComponent(void)
		{
			this->components = (gcnew System::ComponentModel::Container());
			System::Windows::Forms::DataGridViewCellStyle^ dataGridViewCellStyle13 = (gcnew System::Windows::Forms::DataGridViewCellStyle());
			System::Windows::Forms::DataGridViewCellStyle^ dataGridViewCellStyle14 = (gcnew System::Windows::Forms::DataGridViewCellStyle());
			this->zedGraphControl1 = (gcnew ZedGraph::ZedGraphControl());
			this->button1 = (gcnew System::Windows::Forms::Button());
			this->button3 = (gcnew System::Windows::Forms::Button());
			this->dataGridView1 = (gcnew System::Windows::Forms::DataGridView());
			this->i = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->x_i = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->v_i = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->v2_i = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->v_i_v2_i = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->OLP = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->h_i = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->C1 = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->C2 = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->u_i = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->u_i_v_i = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->label1 = (gcnew System::Windows::Forms::Label());
			this->textBox1 = (gcnew System::Windows::Forms::TextBox());
			this->label2 = (gcnew System::Windows::Forms::Label());
			this->textBox2 = (gcnew System::Windows::Forms::TextBox());
			this->label3 = (gcnew System::Windows::Forms::Label());
			this->textBox3 = (gcnew System::Windows::Forms::TextBox());
			this->button2 = (gcnew System::Windows::Forms::Button());
			this->textBox6 = (gcnew System::Windows::Forms::TextBox());
			this->label6 = (gcnew System::Windows::Forms::Label());
			this->textBox7 = (gcnew System::Windows::Forms::TextBox());
			this->label7 = (gcnew System::Windows::Forms::Label());
			this->textBox4 = (gcnew System::Windows::Forms::TextBox());
			this->label4 = (gcnew System::Windows::Forms::Label());
			this->button4 = (gcnew System::Windows::Forms::Button());
			this->label5 = (gcnew System::Windows::Forms::Label());
			this->label8 = (gcnew System::Windows::Forms::Label());
			this->label9 = (gcnew System::Windows::Forms::Label());
			this->button5 = (gcnew System::Windows::Forms::Button());
			this->button6 = (gcnew System::Windows::Forms::Button());
			this->button7 = (gcnew System::Windows::Forms::Button());
			this->label10 = (gcnew System::Windows::Forms::Label());
			this->label11 = (gcnew System::Windows::Forms::Label());
			this->textBox5 = (gcnew System::Windows::Forms::TextBox());
			this->textBox8 = (gcnew System::Windows::Forms::TextBox());
			this->label12 = (gcnew System::Windows::Forms::Label());
			this->label13 = (gcnew System::Windows::Forms::Label());
			this->textBox9 = (gcnew System::Windows::Forms::TextBox());
			this->textBox10 = (gcnew System::Windows::Forms::TextBox());
			this->textBox11 = (gcnew System::Windows::Forms::TextBox());
			this->label14 = (gcnew System::Windows::Forms::Label());
			this->textBox12 = (gcnew System::Windows::Forms::TextBox());
			this->textBox13 = (gcnew System::Windows::Forms::TextBox());
			this->textBox14 = (gcnew System::Windows::Forms::TextBox());
			this->label15 = (gcnew System::Windows::Forms::Label());
			this->label16 = (gcnew System::Windows::Forms::Label());
			this->label17 = (gcnew System::Windows::Forms::Label());
			this->label18 = (gcnew System::Windows::Forms::Label());
			this->label19 = (gcnew System::Windows::Forms::Label());
			this->textBox15 = (gcnew System::Windows::Forms::TextBox());
			this->textBox16 = (gcnew System::Windows::Forms::TextBox());
			this->textBox17 = (gcnew System::Windows::Forms::TextBox());
			this->textBox18 = (gcnew System::Windows::Forms::TextBox());
			this->label20 = (gcnew System::Windows::Forms::Label());
			this->label21 = (gcnew System::Windows::Forms::Label());
			this->label22 = (gcnew System::Windows::Forms::Label());
			this->label23 = (gcnew System::Windows::Forms::Label());
			this->textBox19 = (gcnew System::Windows::Forms::TextBox());
			this->label24 = (gcnew System::Windows::Forms::Label());
			this->label25 = (gcnew System::Windows::Forms::Label());
			this->label26 = (gcnew System::Windows::Forms::Label());
			this->textBox20 = (gcnew System::Windows::Forms::TextBox());
			this->textBox21 = (gcnew System::Windows::Forms::TextBox());
			this->label27 = (gcnew System::Windows::Forms::Label());
			this->textBox22 = (gcnew System::Windows::Forms::TextBox());
			this->label28 = (gcnew System::Windows::Forms::Label());
			this->textBox23 = (gcnew System::Windows::Forms::TextBox());
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->dataGridView1))->BeginInit();
			this->SuspendLayout();
			// 
			// zedGraphControl1
			// 
			this->zedGraphControl1->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->zedGraphControl1->Location = System::Drawing::Point(12, 30);
			this->zedGraphControl1->Margin = System::Windows::Forms::Padding(4, 5, 4, 5);
			this->zedGraphControl1->Name = L"zedGraphControl1";
			this->zedGraphControl1->ScrollGrace = 0;
			this->zedGraphControl1->ScrollMaxX = 0;
			this->zedGraphControl1->ScrollMaxY = 0;
			this->zedGraphControl1->ScrollMaxY2 = 0;
			this->zedGraphControl1->ScrollMinX = 0;
			this->zedGraphControl1->ScrollMinY = 0;
			this->zedGraphControl1->ScrollMinY2 = 0;
			this->zedGraphControl1->Size = System::Drawing::Size(501, 327);
			this->zedGraphControl1->TabIndex = 0;
			this->zedGraphControl1->Load += gcnew System::EventHandler(this, &MyForm::zedGraphControl1_Load);
			// 
			// button1
			// 
			this->button1->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->button1->Location = System::Drawing::Point(600, 395);
			this->button1->Name = L"button1";
			this->button1->Size = System::Drawing::Size(142, 29);
			this->button1->TabIndex = 1;
			this->button1->Text = L"Без контроля";
			this->button1->UseVisualStyleBackColor = true;
			this->button1->Click += gcnew System::EventHandler(this, &MyForm::button1_Click);
			// 
			// button3
			// 
			this->button3->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->button3->Location = System::Drawing::Point(600, 437);
			this->button3->Name = L"button3";
			this->button3->Size = System::Drawing::Size(142, 52);
			this->button3->TabIndex = 3;
			this->button3->Text = L"С оценокой лок погрешности";
			this->button3->UseVisualStyleBackColor = true;
			this->button3->Click += gcnew System::EventHandler(this, &MyForm::button3_Click);
			// 
			// dataGridView1
			// 
			dataGridViewCellStyle13->Alignment = System::Windows::Forms::DataGridViewContentAlignment::MiddleLeft;
			dataGridViewCellStyle13->BackColor = System::Drawing::SystemColors::Control;
			dataGridViewCellStyle13->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 11.25F, System::Drawing::FontStyle::Regular,
				System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(204)));
			dataGridViewCellStyle13->ForeColor = System::Drawing::SystemColors::WindowText;
			dataGridViewCellStyle13->SelectionBackColor = System::Drawing::SystemColors::Highlight;
			dataGridViewCellStyle13->SelectionForeColor = System::Drawing::SystemColors::HighlightText;
			dataGridViewCellStyle13->WrapMode = System::Windows::Forms::DataGridViewTriState::True;
			this->dataGridView1->ColumnHeadersDefaultCellStyle = dataGridViewCellStyle13;
			this->dataGridView1->ColumnHeadersHeightSizeMode = System::Windows::Forms::DataGridViewColumnHeadersHeightSizeMode::AutoSize;
			this->dataGridView1->Columns->AddRange(gcnew cli::array< System::Windows::Forms::DataGridViewColumn^  >(11) {
				this->i, this->x_i,
					this->v_i, this->v2_i, this->v_i_v2_i, this->OLP, this->h_i, this->C1, this->C2, this->u_i, this->u_i_v_i
			});
			dataGridViewCellStyle14->Alignment = System::Windows::Forms::DataGridViewContentAlignment::MiddleLeft;
			dataGridViewCellStyle14->BackColor = System::Drawing::SystemColors::Window;
			dataGridViewCellStyle14->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 11.25F, System::Drawing::FontStyle::Regular,
				System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(204)));
			dataGridViewCellStyle14->ForeColor = System::Drawing::SystemColors::ControlText;
			dataGridViewCellStyle14->SelectionBackColor = System::Drawing::SystemColors::Highlight;
			dataGridViewCellStyle14->SelectionForeColor = System::Drawing::SystemColors::HighlightText;
			dataGridViewCellStyle14->WrapMode = System::Windows::Forms::DataGridViewTriState::False;
			this->dataGridView1->DefaultCellStyle = dataGridViewCellStyle14;
			this->dataGridView1->Location = System::Drawing::Point(519, 30);
			this->dataGridView1->Name = L"dataGridView1";
			this->dataGridView1->RowHeadersVisible = false;
			this->dataGridView1->Size = System::Drawing::Size(1058, 327);
			this->dataGridView1->TabIndex = 2;
			// 
			// i
			// 
			this->i->HeaderText = L"i";
			this->i->Name = L"i";
			this->i->ReadOnly = true;
			this->i->Width = 50;
			// 
			// x_i
			// 
			this->x_i->HeaderText = L"x_i";
			this->x_i->Name = L"x_i";
			this->x_i->ReadOnly = true;
			// 
			// v_i
			// 
			this->v_i->HeaderText = L"v_i";
			this->v_i->Name = L"v_i";
			this->v_i->ReadOnly = true;
			// 
			// v2_i
			// 
			this->v2_i->HeaderText = L"v2_i";
			this->v2_i->Name = L"v2_i";
			this->v2_i->ReadOnly = true;
			// 
			// v_i_v2_i
			// 
			this->v_i_v2_i->HeaderText = L"v_i-v2_i";
			this->v_i_v2_i->Name = L"v_i_v2_i";
			this->v_i_v2_i->ReadOnly = true;
			// 
			// OLP
			// 
			this->OLP->HeaderText = L"ОЛП";
			this->OLP->Name = L"OLP";
			this->OLP->ReadOnly = true;
			this->OLP->Width = 170;
			// 
			// h_i
			// 
			this->h_i->HeaderText = L"h_i";
			this->h_i->Name = L"h_i";
			this->h_i->ReadOnly = true;
			// 
			// C1
			// 
			this->C1->HeaderText = L"C1";
			this->C1->Name = L"C1";
			this->C1->ReadOnly = true;
			// 
			// C2
			// 
			this->C2->HeaderText = L"C2";
			this->C2->Name = L"C2";
			this->C2->ReadOnly = true;
			// 
			// u_i
			// 
			this->u_i->HeaderText = L"u_i";
			this->u_i->Name = L"u_i";
			this->u_i->ReadOnly = true;
			// 
			// u_i_v_i
			// 
			this->u_i_v_i->HeaderText = L"|u_i-v_i|";
			this->u_i_v_i->Name = L"u_i_v_i";
			this->u_i_v_i->ReadOnly = true;
			// 
			// label1
			// 
			this->label1->AutoSize = true;
			this->label1->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label1->Location = System::Drawing::Point(23, 397);
			this->label1->Name = L"label1";
			this->label1->Size = System::Drawing::Size(50, 20);
			this->label1->TabIndex = 3;
			this->label1->Text = L"x_min";
			// 
			// textBox1
			// 
			this->textBox1->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->textBox1->Location = System::Drawing::Point(79, 394);
			this->textBox1->Name = L"textBox1";
			this->textBox1->Size = System::Drawing::Size(48, 26);
			this->textBox1->TabIndex = 4;
			this->textBox1->Text = L"0";
			// 
			// label2
			// 
			this->label2->AutoSize = true;
			this->label2->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label2->Location = System::Drawing::Point(137, 396);
			this->label2->Name = L"label2";
			this->label2->Size = System::Drawing::Size(54, 20);
			this->label2->TabIndex = 5;
			this->label2->Text = L"x_max";
			// 
			// textBox2
			// 
			this->textBox2->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->textBox2->Location = System::Drawing::Point(197, 396);
			this->textBox2->Name = L"textBox2";
			this->textBox2->Size = System::Drawing::Size(49, 26);
			this->textBox2->TabIndex = 6;
			this->textBox2->Text = L"1";
			// 
			// label3
			// 
			this->label3->AutoSize = true;
			this->label3->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label3->Location = System::Drawing::Point(270, 396);
			this->label3->Name = L"label3";
			this->label3->Size = System::Drawing::Size(18, 20);
			this->label3->TabIndex = 7;
			this->label3->Text = L"h";
			// 
			// textBox3
			// 
			this->textBox3->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->textBox3->Location = System::Drawing::Point(294, 395);
			this->textBox3->Name = L"textBox3";
			this->textBox3->Size = System::Drawing::Size(61, 26);
			this->textBox3->TabIndex = 8;
			this->textBox3->Text = L"0,1";
			// 
			// button2
			// 
			this->button2->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->button2->Location = System::Drawing::Point(1348, 440);
			this->button2->Name = L"button2";
			this->button2->Size = System::Drawing::Size(141, 59);
			this->button2->TabIndex = 9;
			this->button2->Text = L"Фазовый портрет";
			this->button2->UseVisualStyleBackColor = true;
			this->button2->Click += gcnew System::EventHandler(this, &MyForm::button2_Click);
			// 
			// textBox6
			// 
			this->textBox6->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->textBox6->Location = System::Drawing::Point(513, 397);
			this->textBox6->Name = L"textBox6";
			this->textBox6->Size = System::Drawing::Size(61, 26);
			this->textBox6->TabIndex = 14;
			this->textBox6->Text = L"1";
			// 
			// label6
			// 
			this->label6->AutoSize = true;
			this->label6->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label6->Location = System::Drawing::Point(480, 400);
			this->label6->Name = L"label6";
			this->label6->Size = System::Drawing::Size(27, 20);
			this->label6->TabIndex = 15;
			this->label6->Text = L"u0";
			// 
			// textBox7
			// 
			this->textBox7->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->textBox7->Location = System::Drawing::Point(385, 396);
			this->textBox7->Name = L"textBox7";
			this->textBox7->Size = System::Drawing::Size(89, 26);
			this->textBox7->TabIndex = 17;
			this->textBox7->Text = L"0,00001";
			// 
			// label7
			// 
			this->label7->AutoSize = true;
			this->label7->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label7->Location = System::Drawing::Point(361, 398);
			this->label7->Name = L"label7";
			this->label7->Size = System::Drawing::Size(18, 20);
			this->label7->TabIndex = 18;
			this->label7->Text = L"e";
			// 
			// textBox4
			// 
			this->textBox4->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->textBox4->Location = System::Drawing::Point(1428, 400);
			this->textBox4->Name = L"textBox4";
			this->textBox4->Size = System::Drawing::Size(61, 26);
			this->textBox4->TabIndex = 19;
			this->textBox4->Text = L"0,01";
			// 
			// label4
			// 
			this->label4->AutoSize = true;
			this->label4->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label4->Location = System::Drawing::Point(1350, 399);
			this->label4->Name = L"label4";
			this->label4->Size = System::Drawing::Size(72, 20);
			this->label4->TabIndex = 20;
			this->label4->Text = L"Граница";
			// 
			// button4
			// 
			this->button4->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->button4->Location = System::Drawing::Point(748, 437);
			this->button4->Name = L"button4";
			this->button4->Size = System::Drawing::Size(142, 52);
			this->button4->TabIndex = 21;
			this->button4->Text = L"С оценокой лок погрешности";
			this->button4->UseVisualStyleBackColor = true;
			this->button4->Click += gcnew System::EventHandler(this, &MyForm::button4_Click);
			// 
			// label5
			// 
			this->label5->AutoSize = true;
			this->label5->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label5->Location = System::Drawing::Point(603, 372);
			this->label5->Name = L"label5";
			this->label5->Size = System::Drawing::Size(139, 20);
			this->label5->TabIndex = 22;
			this->label5->Text = L"Тестовая задача";
			// 
			// label8
			// 
			this->label8->AutoSize = true;
			this->label8->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label8->Location = System::Drawing::Point(785, 372);
			this->label8->Name = L"label8";
			this->label8->Size = System::Drawing::Size(80, 20);
			this->label8->TabIndex = 23;
			this->label8->Text = L"Задача 1";
			// 
			// label9
			// 
			this->label9->AutoSize = true;
			this->label9->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label9->Location = System::Drawing::Point(921, 372);
			this->label9->Name = L"label9";
			this->label9->Size = System::Drawing::Size(80, 20);
			this->label9->TabIndex = 24;
			this->label9->Text = L"Задача 2";
			// 
			// button5
			// 
			this->button5->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->button5->Location = System::Drawing::Point(748, 395);
			this->button5->Name = L"button5";
			this->button5->Size = System::Drawing::Size(142, 29);
			this->button5->TabIndex = 25;
			this->button5->Text = L"Без контроля";
			this->button5->UseVisualStyleBackColor = true;
			this->button5->Click += gcnew System::EventHandler(this, &MyForm::button5_Click);
			// 
			// button6
			// 
			this->button6->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->button6->Location = System::Drawing::Point(896, 396);
			this->button6->Name = L"button6";
			this->button6->Size = System::Drawing::Size(142, 29);
			this->button6->TabIndex = 26;
			this->button6->Text = L"Без контроля";
			this->button6->UseVisualStyleBackColor = true;
			this->button6->Click += gcnew System::EventHandler(this, &MyForm::button6_Click);
			// 
			// button7
			// 
			this->button7->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->button7->Location = System::Drawing::Point(896, 437);
			this->button7->Name = L"button7";
			this->button7->Size = System::Drawing::Size(142, 52);
			this->button7->TabIndex = 27;
			this->button7->Text = L"С оценокой лок погрешности";
			this->button7->UseVisualStyleBackColor = true;
			this->button7->Click += gcnew System::EventHandler(this, &MyForm::button7_Click);
			// 
			// label10
			// 
			this->label10->AutoSize = true;
			this->label10->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label10->Location = System::Drawing::Point(1060, 399);
			this->label10->Name = L"label10";
			this->label10->Size = System::Drawing::Size(45, 20);
			this->label10->TabIndex = 28;
			this->label10->Text = L"u0_1";
			// 
			// label11
			// 
			this->label11->AutoSize = true;
			this->label11->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label11->Location = System::Drawing::Point(1060, 440);
			this->label11->Name = L"label11";
			this->label11->Size = System::Drawing::Size(45, 20);
			this->label11->TabIndex = 29;
			this->label11->Text = L"u0_2";
			// 
			// textBox5
			// 
			this->textBox5->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->textBox5->Location = System::Drawing::Point(1111, 395);
			this->textBox5->Name = L"textBox5";
			this->textBox5->Size = System::Drawing::Size(61, 26);
			this->textBox5->TabIndex = 30;
			this->textBox5->Text = L"1";
			// 
			// textBox8
			// 
			this->textBox8->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->textBox8->Location = System::Drawing::Point(1111, 437);
			this->textBox8->Name = L"textBox8";
			this->textBox8->Size = System::Drawing::Size(61, 26);
			this->textBox8->TabIndex = 31;
			this->textBox8->Text = L"1";
			// 
			// label12
			// 
			this->label12->AutoSize = true;
			this->label12->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label12->Location = System::Drawing::Point(1191, 400);
			this->label12->Name = L"label12";
			this->label12->Size = System::Drawing::Size(18, 20);
			this->label12->TabIndex = 32;
			this->label12->Text = L"a";
			// 
			// label13
			// 
			this->label13->AutoSize = true;
			this->label13->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label13->Location = System::Drawing::Point(1191, 443);
			this->label13->Name = L"label13";
			this->label13->Size = System::Drawing::Size(18, 20);
			this->label13->TabIndex = 33;
			this->label13->Text = L"b";
			// 
			// textBox9
			// 
			this->textBox9->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->textBox9->Location = System::Drawing::Point(1215, 395);
			this->textBox9->Name = L"textBox9";
			this->textBox9->Size = System::Drawing::Size(61, 26);
			this->textBox9->TabIndex = 34;
			this->textBox9->Text = L"1";
			// 
			// textBox10
			// 
			this->textBox10->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->textBox10->Location = System::Drawing::Point(1215, 437);
			this->textBox10->Name = L"textBox10";
			this->textBox10->Size = System::Drawing::Size(61, 26);
			this->textBox10->TabIndex = 35;
			this->textBox10->Text = L"1";
			// 
			// textBox11
			// 
			this->textBox11->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->textBox11->Location = System::Drawing::Point(1428, 366);
			this->textBox11->Name = L"textBox11";
			this->textBox11->Size = System::Drawing::Size(61, 26);
			this->textBox11->TabIndex = 36;
			this->textBox11->Text = L"1000";
			// 
			// label14
			// 
			this->label14->AutoSize = true;
			this->label14->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label14->Location = System::Drawing::Point(1286, 372);
			this->label14->Name = L"label14";
			this->label14->Size = System::Drawing::Size(136, 20);
			this->label14->TabIndex = 37;
			this->label14->Text = L"Максимум шагов";
			// 
			// textBox12
			// 
			this->textBox12->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->textBox12->Location = System::Drawing::Point(64, 466);
			this->textBox12->Name = L"textBox12";
			this->textBox12->Size = System::Drawing::Size(89, 26);
			this->textBox12->TabIndex = 38;
			// 
			// textBox13
			// 
			this->textBox13->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->textBox13->Location = System::Drawing::Point(93, 497);
			this->textBox13->Name = L"textBox13";
			this->textBox13->Size = System::Drawing::Size(183, 26);
			this->textBox13->TabIndex = 39;
			// 
			// textBox14
			// 
			this->textBox14->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->textBox14->Location = System::Drawing::Point(113, 532);
			this->textBox14->Name = L"textBox14";
			this->textBox14->Size = System::Drawing::Size(163, 26);
			this->textBox14->TabIndex = 40;
			// 
			// label15
			// 
			this->label15->AutoSize = true;
			this->label15->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label15->Location = System::Drawing::Point(23, 469);
			this->label15->Name = L"label15";
			this->label15->Size = System::Drawing::Size(35, 20);
			this->label15->TabIndex = 41;
			this->label15->Text = L"n = ";
			// 
			// label16
			// 
			this->label16->AutoSize = true;
			this->label16->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label16->Location = System::Drawing::Point(23, 503);
			this->label16->Name = L"label16";
			this->label16->Size = System::Drawing::Size(56, 20);
			this->label16->TabIndex = 42;
			this->label16->Text = L"b - x_n";
			// 
			// label17
			// 
			this->label17->AutoSize = true;
			this->label17->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label17->Location = System::Drawing::Point(23, 532);
			this->label17->Name = L"label17";
			this->label17->Size = System::Drawing::Size(84, 20);
			this->label17->TabIndex = 43;
			this->label17->Text = L"max|ОЛП|";
			// 
			// label18
			// 
			this->label18->AutoSize = true;
			this->label18->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label18->Location = System::Drawing::Point(128, 437);
			this->label18->Name = L"label18";
			this->label18->Size = System::Drawing::Size(238, 20);
			this->label18->TabIndex = 44;
			this->label18->Text = L"Выходные данные программы";
			// 
			// label19
			// 
			this->label19->AutoSize = true;
			this->label19->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label19->Location = System::Drawing::Point(128, 370);
			this->label19->Name = L"label19";
			this->label19->Size = System::Drawing::Size(227, 20);
			this->label19->TabIndex = 45;
			this->label19->Text = L"Входные данные программы";
			// 
			// textBox15
			// 
			this->textBox15->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->textBox15->Location = System::Drawing::Point(347, 502);
			this->textBox15->Name = L"textBox15";
			this->textBox15->Size = System::Drawing::Size(122, 26);
			this->textBox15->TabIndex = 46;
			// 
			// textBox16
			// 
			this->textBox16->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->textBox16->Location = System::Drawing::Point(550, 503);
			this->textBox16->Name = L"textBox16";
			this->textBox16->Size = System::Drawing::Size(90, 26);
			this->textBox16->TabIndex = 47;
			// 
			// textBox17
			// 
			this->textBox17->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->textBox17->Location = System::Drawing::Point(347, 531);
			this->textBox17->Name = L"textBox17";
			this->textBox17->Size = System::Drawing::Size(122, 26);
			this->textBox17->TabIndex = 48;
			// 
			// textBox18
			// 
			this->textBox18->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->textBox18->Location = System::Drawing::Point(550, 535);
			this->textBox18->Name = L"textBox18";
			this->textBox18->Size = System::Drawing::Size(90, 26);
			this->textBox18->TabIndex = 49;
			// 
			// label20
			// 
			this->label20->AutoSize = true;
			this->label20->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label20->Location = System::Drawing::Point(290, 506);
			this->label20->Name = L"label20";
			this->label20->Size = System::Drawing::Size(51, 20);
			this->label20->TabIndex = 50;
			this->label20->Text = L"max h";
			// 
			// label21
			// 
			this->label21->AutoSize = true;
			this->label21->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label21->Location = System::Drawing::Point(294, 535);
			this->label21->Name = L"label21";
			this->label21->Size = System::Drawing::Size(47, 20);
			this->label21->TabIndex = 51;
			this->label21->Text = L"min h";
			// 
			// label22
			// 
			this->label22->AutoSize = true;
			this->label22->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label22->Location = System::Drawing::Point(480, 506);
			this->label22->Name = L"label22";
			this->label22->Size = System::Drawing::Size(64, 20);
			this->label22->TabIndex = 52;
			this->label22->Text = L"при х = ";
			// 
			// label23
			// 
			this->label23->AutoSize = true;
			this->label23->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label23->Location = System::Drawing::Point(23, 567);
			this->label23->Name = L"label23";
			this->label23->Size = System::Drawing::Size(64, 20);
			this->label23->TabIndex = 53;
			this->label23->Text = L"при х = ";
			// 
			// textBox19
			// 
			this->textBox19->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->textBox19->Location = System::Drawing::Point(79, 565);
			this->textBox19->Name = L"textBox19";
			this->textBox19->Size = System::Drawing::Size(89, 26);
			this->textBox19->TabIndex = 54;
			// 
			// label24
			// 
			this->label24->AutoSize = true;
			this->label24->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label24->Location = System::Drawing::Point(480, 537);
			this->label24->Name = L"label24";
			this->label24->Size = System::Drawing::Size(64, 20);
			this->label24->TabIndex = 55;
			this->label24->Text = L"при х = ";
			// 
			// label25
			// 
			this->label25->AutoSize = true;
			this->label25->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label25->Location = System::Drawing::Point(673, 506);
			this->label25->Name = L"label25";
			this->label25->Size = System::Drawing::Size(192, 20);
			this->label25->TabIndex = 56;
			this->label25->Text = L"Число увеличений шага";
			// 
			// label26
			// 
			this->label26->AutoSize = true;
			this->label26->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label26->Location = System::Drawing::Point(673, 541);
			this->label26->Name = L"label26";
			this->label26->Size = System::Drawing::Size(197, 20);
			this->label26->TabIndex = 57;
			this->label26->Text = L"Число уменьшений шага";
			// 
			// textBox20
			// 
			this->textBox20->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->textBox20->Location = System::Drawing::Point(871, 505);
			this->textBox20->Name = L"textBox20";
			this->textBox20->Size = System::Drawing::Size(67, 26);
			this->textBox20->TabIndex = 58;
			// 
			// textBox21
			// 
			this->textBox21->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->textBox21->Location = System::Drawing::Point(871, 541);
			this->textBox21->Name = L"textBox21";
			this->textBox21->Size = System::Drawing::Size(67, 26);
			this->textBox21->TabIndex = 59;
			// 
			// label27
			// 
			this->label27->AutoSize = true;
			this->label27->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label27->Location = System::Drawing::Point(215, 571);
			this->label27->Name = L"label27";
			this->label27->Size = System::Drawing::Size(93, 20);
			this->label27->TabIndex = 60;
			this->label27->Text = L"max|u_i-v_i|";
			// 
			// textBox22
			// 
			this->textBox22->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->textBox22->Location = System::Drawing::Point(314, 568);
			this->textBox22->Name = L"textBox22";
			this->textBox22->Size = System::Drawing::Size(155, 26);
			this->textBox22->TabIndex = 61;
			// 
			// label28
			// 
			this->label28->AutoSize = true;
			this->label28->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label28->Location = System::Drawing::Point(480, 571);
			this->label28->Name = L"label28";
			this->label28->Size = System::Drawing::Size(64, 20);
			this->label28->TabIndex = 62;
			this->label28->Text = L"при х = ";
			// 
			// textBox23
			// 
			this->textBox23->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->textBox23->Location = System::Drawing::Point(550, 568);
			this->textBox23->Name = L"textBox23";
			this->textBox23->Size = System::Drawing::Size(90, 26);
			this->textBox23->TabIndex = 63;
			// 
			// MyForm
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->ClientSize = System::Drawing::Size(1576, 602);
			this->Controls->Add(this->textBox23);
			this->Controls->Add(this->label28);
			this->Controls->Add(this->textBox22);
			this->Controls->Add(this->label27);
			this->Controls->Add(this->textBox21);
			this->Controls->Add(this->textBox20);
			this->Controls->Add(this->label26);
			this->Controls->Add(this->label25);
			this->Controls->Add(this->label24);
			this->Controls->Add(this->textBox19);
			this->Controls->Add(this->label23);
			this->Controls->Add(this->label22);
			this->Controls->Add(this->label21);
			this->Controls->Add(this->label20);
			this->Controls->Add(this->textBox18);
			this->Controls->Add(this->textBox17);
			this->Controls->Add(this->textBox16);
			this->Controls->Add(this->textBox15);
			this->Controls->Add(this->label19);
			this->Controls->Add(this->label18);
			this->Controls->Add(this->label17);
			this->Controls->Add(this->label16);
			this->Controls->Add(this->label15);
			this->Controls->Add(this->textBox14);
			this->Controls->Add(this->textBox13);
			this->Controls->Add(this->textBox12);
			this->Controls->Add(this->label14);
			this->Controls->Add(this->textBox11);
			this->Controls->Add(this->textBox10);
			this->Controls->Add(this->textBox9);
			this->Controls->Add(this->label13);
			this->Controls->Add(this->label12);
			this->Controls->Add(this->textBox8);
			this->Controls->Add(this->textBox5);
			this->Controls->Add(this->label11);
			this->Controls->Add(this->label10);
			this->Controls->Add(this->button7);
			this->Controls->Add(this->button6);
			this->Controls->Add(this->button5);
			this->Controls->Add(this->label9);
			this->Controls->Add(this->label8);
			this->Controls->Add(this->label5);
			this->Controls->Add(this->button4);
			this->Controls->Add(this->label4);
			this->Controls->Add(this->textBox4);
			this->Controls->Add(this->label7);
			this->Controls->Add(this->textBox7);
			this->Controls->Add(this->label6);
			this->Controls->Add(this->textBox6);
			this->Controls->Add(this->button2);
			this->Controls->Add(this->textBox3);
			this->Controls->Add(this->label3);
			this->Controls->Add(this->textBox2);
			this->Controls->Add(this->label2);
			this->Controls->Add(this->textBox1);
			this->Controls->Add(this->label1);
			this->Controls->Add(this->dataGridView1);
			this->Controls->Add(this->button1);
			this->Controls->Add(this->button3);
			this->Controls->Add(this->zedGraphControl1);
			this->Name = L"MyForm";
			this->Text = L"MyForm";
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->dataGridView1))->EndInit();
			this->ResumeLayout(false);
			this->PerformLayout();

		}
#pragma endregion
// Тестовая задача без контроля локальной погрешности
	private: System::Void button1_Click(System::Object^ sender, System::EventArgs^ e) {

		GraphPane^ panel = zedGraphControl1->GraphPane;
		panel->CurveList->Clear();
		PointPairList^ true_list = gcnew ZedGraph::PointPairList();
		PointPairList^ test_list = gcnew ZedGraph::PointPairList();

		// Интервал, где есть данные
		double xmin = Convert::ToDouble(textBox1->Text);
		double xmax = Convert::ToDouble(textBox2->Text);
		double border = Convert::ToDouble(textBox4->Text);
		double h = Convert::ToDouble(textBox3->Text);
		double u0 = Convert::ToDouble(textBox6->Text);
		std::size_t Max_steps = Convert::ToInt64(textBox11->Text);

		textBox12->Clear();
		textBox13->Clear();
		textBox14->Clear();
		textBox15->Clear();
		textBox16->Clear();
		textBox17->Clear();
		textBox18->Clear();
		textBox19->Clear();
		textBox20->Clear();
		textBox21->Clear();
		textBox22->Clear();
		textBox23->Clear();

		double xmin_limit = xmin - 0.1;
		double xmax_limit = xmax + 0.1;

		// Список точек
		int i = 1;
		dataGridView1->Rows->Clear();
		double v = u0;
		double v_last = u0;
		double v_2h;
		double v_ist;

		dataGridView1->Rows->Add();
		dataGridView1->Rows[0]->Cells[0]->Value = 0.0;
		dataGridView1->Rows[0]->Cells[1]->Value = xmin;
		dataGridView1->Rows[0]->Cells[2]->Value = floor(v * 10000) / 10000;
		dataGridView1->Rows[0]->Cells[9]->Value = u0;
		dataGridView1->Rows[0]->Cells[10]->Value = abs(u0 - v);

		test_list->Add(xmin, v);
		true_list->Add(xmin, true_trajectory(xmin, u0));

		double x = xmin + h;
		double max_OLP = 0.0;
		double max_OLP_x = 0.0;
		double min_h = h;
		double max_h = 0;
		double min_h_x = 0, max_h_x = 0;
		double max_v_ist_v = 0; double max_v_ist_v_x = 0;

		for (; (x <= xmax) && (i < Max_steps); x += h)
		{
			//Добавление на график
			true_list->Add(x, true_trajectory(x, u0));
			//f2_list->Add(x, f2(x));
			v_last = v;
			v = Runge_Kytta_4(test_function, h, x, v).second;
			test_list->Add(x, v);

			//Печать в таблицу
			v_2h = Runge_Kytta_4(test_function, h / 2.0, x + h / 2.0, Runge_Kytta_4(test_function, h / 2.0, x, v_last).second).second;
			v_ist = floor(true_trajectory(x, u0) * 10000) / 10000;
			dataGridView1->Rows->Add();
			dataGridView1->Rows[i]->Cells[0]->Value = i;
			dataGridView1->Rows[i]->Cells[1]->Value = x;
			dataGridView1->Rows[i]->Cells[2]->Value = floor(v * 10000) / 10000;
			dataGridView1->Rows[i]->Cells[3]->Value = floor(v_2h * 10000) / 10000;
			dataGridView1->Rows[i]->Cells[4]->Value = abs(v - v_2h);
			dataGridView1->Rows[i]->Cells[5]->Value = S(v, v_2h) * pow(2, p); if (S(v, v_2h) * pow(2, p) > max_OLP) { max_OLP = S(v, v_2h) * pow(2, p); max_OLP_x = x; }
			dataGridView1->Rows[i]->Cells[6]->Value = h; if (h > max_h) { max_h = h; max_h_x = x; } if (h < min_h) { min_h = h; min_h_x = x; }
			dataGridView1->Rows[i]->Cells[7]->Value = 0.0;
			dataGridView1->Rows[i]->Cells[8]->Value = 0.0;
			dataGridView1->Rows[i]->Cells[9]->Value = v_ist;
			dataGridView1->Rows[i]->Cells[10]->Value = abs(v_ist - v); if (abs(v_ist - v) > max_v_ist_v) { max_v_ist_v = abs(v_ist - v); max_v_ist_v_x = x; }
			i++;
		}
		LineItem^ Curve1 = panel->AddCurve("true_trajectory(x)", true_list, Color::Red, SymbolType::Plus);
		LineItem^ Curve2 = panel->AddCurve("test_function(x)", test_list, Color::Blue, SymbolType::None);

		textBox12->AppendText(Convert::ToString(i));
		textBox13->AppendText(Convert::ToString(xmax - x+h));
		textBox14->AppendText(Convert::ToString(max_OLP));
		textBox15->AppendText(Convert::ToString(max_h));
		textBox16->AppendText(Convert::ToString(max_h_x));
		textBox17->AppendText(Convert::ToString(min_h));
		textBox18->AppendText(Convert::ToString(min_h_x));
		textBox19->AppendText(Convert::ToString(max_OLP_x));
		textBox20->AppendText(Convert::ToString(0));
		textBox21->AppendText(Convert::ToString(0));
		textBox22->AppendText(Convert::ToString(max_v_ist_v));
		textBox23->AppendText(Convert::ToString(max_v_ist_v_x));


		// Устанавливаем интересующий нас интервал по оси X
		panel->XAxis->Scale->Min = xmin_limit;
		panel->XAxis->Scale->Max = xmax_limit;
		/*
				// Устанавливаем интересующий нас интервал по оси Y
				panel->YAxis->Scale->Min = ymin_limit;
				panel->YAxis->Scale->Max = ymax_limit;
		*/
		// Вызываем метод AxisChange (), чтобы обновить данные об осях. 
		// В противном случае на рисунке будет показана только часть графика, 
		// которая умещается в интервалы по осям, установленные по умолчанию
		zedGraphControl1->AxisChange();
		// Обновляем график
		zedGraphControl1->Invalidate();

	}
// Тестовая задача с контролем локальной погрешности
	private: System::Void button3_Click(System::Object^ sender, System::EventArgs^ e) {

		GraphPane^ panel = zedGraphControl1->GraphPane;
		panel->CurveList->Clear();
		PointPairList^ true_list = gcnew ZedGraph::PointPairList();
		PointPairList^ test_list = gcnew ZedGraph::PointPairList();

		// Интервал, где есть данные
		double xmin = Convert::ToDouble(textBox1->Text);
		double xmax = Convert::ToDouble(textBox2->Text);
		double e_0 = Convert::ToDouble(textBox7->Text);
		double border = Convert::ToDouble(textBox4->Text);
		std::size_t Max_steps = Convert::ToInt64(textBox11->Text);

		textBox12->Clear();
		textBox13->Clear();
		textBox14->Clear();
		textBox15->Clear();
		textBox16->Clear();
		textBox17->Clear();
		textBox18->Clear();
		textBox19->Clear();
		textBox20->Clear();
		textBox21->Clear();
		textBox22->Clear();
		textBox23->Clear();

		double h = Convert::ToDouble(textBox3->Text);
		double u0 = Convert::ToDouble(textBox6->Text);

		double xmin_limit = xmin - 0.1;
		double xmax_limit = xmax + 0.1;
		// Список точек
		int i = 1;
		dataGridView1->Rows->Clear();
		double v = u0;
		double v_last = u0;
		double last_h = 0;
		double v_2h;
		double v_ist;
		std::size_t C1 = 0, C2 = 0;
		std::size_t C1_amount = 0, C2_amount = 0;

		//test_list->Add(xmin, u0);
		std::vector<double> new_point, last_point;
		new_point = { xmin, v, h };

		test_list->Add(xmin, u0);
		true_list->Add(xmin, true_trajectory(xmin, u0));

		dataGridView1->Rows->Add();
		dataGridView1->Rows[0]->Cells[0]->Value = 0.0;
		dataGridView1->Rows[0]->Cells[1]->Value = xmin;
		dataGridView1->Rows[0]->Cells[2]->Value = floor(v * 10000) / 10000;
		dataGridView1->Rows[0]->Cells[9]->Value = u0;
		dataGridView1->Rows[0]->Cells[10]->Value = abs(u0 - v);

		double x = xmin;
		double max_OLP = 0.0;
		double max_OLP_x = 0.0;
		double min_h = h;
		double max_h = 0;
		double min_h_x = 0, max_h_x = 0;
		double max_v_ist_v = 0; double max_v_ist_v_x = 0;

		for (; (x < xmax - border) && (i < Max_steps); ) {

			//Добавление на график
			last_h = h;
			last_point = new_point;

			new_point = RK_4_OLP(test_function, x, v, h, e_0);

			x = new_point[0];

			if (x > xmax) {
				new_point = RK_4_OLP(test_function, last_point[0], v, xmax - last_point[0], e_0);
				x = new_point[0];
				if (xmax - last_point[0] - h > last_h) {
					C1 = 0; C2 = 0;
				}
			}
			else {
				if (new_point[4] < 0)
					C1 = new_point[4] * -1;
				if (new_point[4] > 0) {
					C2 = 1;
					if (x + h > xmax)
						C2 = 0;
				}
			}

			x = new_point[0];
			v = new_point[1];
			v_2h = new_point[2];
			h = new_point[3];

			test_list->Add(x, v);
			true_list->Add(x, true_trajectory(x, u0));
			//Печать в таблицу
			v_ist = true_trajectory(x, u0);
			dataGridView1->Rows->Add();
			dataGridView1->Rows[i]->Cells[0]->Value = i;
			dataGridView1->Rows[i]->Cells[1]->Value = x;
			dataGridView1->Rows[i]->Cells[2]->Value = v;
			dataGridView1->Rows[i]->Cells[3]->Value = v_2h;
			dataGridView1->Rows[i]->Cells[4]->Value = abs(v - v_2h);
			dataGridView1->Rows[i]->Cells[5]->Value = abs(v - v_2h) * pow(2, p); if (abs(v - v_2h) * pow(2, p) > max_OLP) { max_OLP = abs(v - v_2h) * pow(2, p); max_OLP_x = x; }
			dataGridView1->Rows[i]->Cells[6]->Value = x - last_point[0]; if (x - last_point[0] > max_h) { max_h = x - last_point[0]; max_h_x = x; } if (x - last_point[0] < min_h) { min_h = x - last_point[0]; min_h_x = x; }
			dataGridView1->Rows[i]->Cells[7]->Value = C1;
			dataGridView1->Rows[i]->Cells[8]->Value = C2;
			dataGridView1->Rows[i]->Cells[9]->Value = v_ist;
			dataGridView1->Rows[i]->Cells[10]->Value = abs(v_ist - v); if (abs(v_ist - v) > max_v_ist_v) { max_v_ist_v = abs(v_ist - v); max_v_ist_v_x = x; }
			i++;
			C1_amount += C1;
			C2_amount += C2;
			C1 = 0;
			C2 = 0;
		}
		LineItem^ Curve1 = panel->AddCurve("true_trajectory(x)", true_list, Color::Red, SymbolType::Plus);
		LineItem^ Curve2 = panel->AddCurve("test_function(x)", test_list, Color::Blue, SymbolType::None);

		textBox12->AppendText(Convert::ToString(i));
		textBox13->AppendText(Convert::ToString(xmax - x));
		textBox14->AppendText(Convert::ToString(max_OLP));
		textBox15->AppendText(Convert::ToString(max_h));
		textBox16->AppendText(Convert::ToString(max_h_x));
		textBox17->AppendText(Convert::ToString(min_h));
		textBox18->AppendText(Convert::ToString(min_h_x));
		textBox19->AppendText(Convert::ToString(max_OLP_x));
		textBox20->AppendText(Convert::ToString(C2_amount));
		textBox21->AppendText(Convert::ToString(C1_amount));
		textBox22->AppendText(Convert::ToString(max_v_ist_v));
		textBox23->AppendText(Convert::ToString(max_v_ist_v_x));


		// Устанавливаем интересующий нас интервал по оси X
		panel->XAxis->Scale->Min = xmin_limit;
		panel->XAxis->Scale->Max = xmax_limit;
		/*
				// Устанавливаем интересующий нас интервал по оси Y
				panel->YAxis->Scale->Min = ymin_limit;
				panel->YAxis->Scale->Max = ymax_limit;
		*/
		// Вызываем метод AxisChange (), чтобы обновить данные об осях. 
		// В противном случае на рисунке будет показана только часть графика, 
		// которая умещается в интервалы по осям, установленные по умолчанию
		zedGraphControl1->AxisChange();
		// Обновляем график
		zedGraphControl1->Invalidate();
	}
// Фазовый портрет для системы дифференциальных уравнений задачи 2
	private: System::Void button2_Click(System::Object^ sender, System::EventArgs^ e) {

		GraphPane^ panel = zedGraphControl1->GraphPane;
		panel->CurveList->Clear();
		PointPairList^ phase_list = gcnew ZedGraph::PointPairList();

		textBox12->Clear();
		textBox13->Clear();
		textBox14->Clear();
		textBox15->Clear();
		textBox16->Clear();
		textBox17->Clear();
		textBox18->Clear();
		textBox19->Clear();
		textBox20->Clear();
		textBox21->Clear();
		textBox22->Clear();
		textBox23->Clear();


		// Интервал, где есть данные
		double xmin = Convert::ToDouble(textBox1->Text);
		double xmax = Convert::ToDouble(textBox2->Text);
		double e_0 = Convert::ToDouble(textBox7->Text);
		double border = Convert::ToDouble(textBox4->Text);
		double u0_1 = Convert::ToDouble(textBox5->Text);
		double u0_2 = Convert::ToDouble(textBox8->Text);
		double a = Convert::ToDouble(textBox9->Text);
		double b = Convert::ToDouble(textBox10->Text);
		std::size_t Max_steps = Convert::ToInt64(textBox11->Text);

		double h = Convert::ToDouble(textBox3->Text);

		double xmin_limit = xmin - 0.1;
		double xmax_limit = xmax + 0.1;

		// Список точек
		int i = 1;
		dataGridView1->Rows->Clear();
		double v_1 = u0_1;
		double v_2 = u0_2;
		double v_last_1 = u0_1;
		double v_last_2 = u0_2;
		double v_2h_1;
		double v_2h_2;
		std::size_t C1 = 0, C2 = 0, last_h = 0;

		//test_list->Add(xmin, u0);
		std::vector<double> new_point, last_point;
		new_point = { xmin, v_1, v_2, h };

		phase_list->Add(v_1, v_2);

		double x = xmin;
		double max_OLP = 0.0;
		double min_h = h;
		double max_h = 0;
		double min_h_x = 0, max_h_x = 0;

		for (; (x <= xmax - border) && (i < Max_steps); ) {
			//Добавление на график
			last_h = h;
			last_point = new_point;

			new_point = RK_4_OLP_for_system(function_2, x, v_1, v_2, h, e_0, a, b);

			x = new_point[0];
			v_1 = new_point[1];
			v_2 = new_point[2];
			v_2h_1 = new_point[3];
			v_2h_2 = new_point[4];
			h = new_point[5];

			phase_list->Add(v_1, v_2);
			i++;
		}
		LineItem^ Curve1 = panel->AddCurve("phase_portrait", phase_list, Color::Green, SymbolType::None);

		textBox12->AppendText(Convert::ToString(i));

		// Устанавливаем интересующий нас интервал по оси X
		panel->XAxis->Scale->Min = xmin_limit;
		panel->XAxis->Scale->Max = xmax_limit;
		/*
				// Устанавливаем интересующий нас интервал по оси Y
				panel->YAxis->Scale->Min = ymin_limit;
				panel->YAxis->Scale->Max = ymax_limit;
		*/
		// Вызываем метод AxisChange (), чтобы обновить данные об осях. 
		// В противном случае на рисунке будет показана только часть графика, 
		// которая умещается в интервалы по осям, установленные по умолчанию
		zedGraphControl1->AxisChange();
		// Обновляем график
		zedGraphControl1->Invalidate();

	}
// Задача 1 с контролем локальной погрешности
	private: System::Void button4_Click(System::Object^ sender, System::EventArgs^ e) {

		GraphPane^ panel = zedGraphControl1->GraphPane;
		panel->CurveList->Clear();
		PointPairList^ function_list = gcnew ZedGraph::PointPairList();

		// Интервал, где есть данные
		double xmin = Convert::ToDouble(textBox1->Text);
		double xmax = Convert::ToDouble(textBox2->Text);
		double e_0 = Convert::ToDouble(textBox7->Text);
		double border = Convert::ToDouble(textBox4->Text);
		std::size_t Max_steps = Convert::ToInt64(textBox11->Text);

		textBox12->Clear();
		textBox13->Clear();
		textBox14->Clear();
		textBox15->Clear();
		textBox16->Clear();
		textBox17->Clear();
		textBox18->Clear();
		textBox19->Clear();
		textBox20->Clear();
		textBox21->Clear();
		textBox22->Clear();
		textBox23->Clear();

		double h = Convert::ToDouble(textBox3->Text);
		double u0 = Convert::ToDouble(textBox6->Text);

		double xmin_limit = xmin - 0.1;
		double xmax_limit = xmax + 0.1;

		// Список точек
		int i = 1;
		dataGridView1->Rows->Clear();
		double v = u0;
		double v_last = u0;
		double v_2h;
		std::size_t C1 = 0, C2 = 0, last_h = 0;
		std::size_t C1_amount = 0, C2_amount = 0;

		//test_list->Add(xmin, u0);
		std::vector<double> new_point, last_point;
		new_point = { xmin, v, h };

		function_list->Add(xmin, u0);

		dataGridView1->Rows->Add();
		dataGridView1->Rows[0]->Cells[0]->Value = 0.0;
		dataGridView1->Rows[0]->Cells[1]->Value = xmin;
		dataGridView1->Rows[0]->Cells[2]->Value = floor(v * 10000) / 10000;
		//dataGridView1->Rows[0]->Cells[9]->Value = u0;
		//dataGridView1->Rows[0]->Cells[10]->Value = abs(u0 - v);

		double x = xmin;
		double max_OLP = 0;
		double max_OLP_x = 0;
		double min_h = h;
		double max_h = 0;
		double min_h_x = 0, max_h_x = 0;

		for (; (x < xmax - border) && (i < Max_steps); ) {
			//Добавление на график
			last_h = h;
			last_point = new_point;

			new_point = RK_4_OLP(function_1, x, v, h, e_0);

			x = new_point[0];


			if (x > xmax) {
				new_point = RK_4_OLP(function_1, last_point[0], v, xmax - last_point[0], e_0);
				x = new_point[0];
				if (xmax - last_point[0] - h > last_h) {
					C1 = 0; C2 = 0;
				}
			}
			else {
				if (new_point[4] < 0)
					C1 = new_point[4] * -1;
				if (new_point[4] > 0) {
					C2 = 1;
					if (x + h > xmax)
						C2 = 0;
				}
			}

			x = new_point[0];
			v = new_point[1];
			v_2h = new_point[2];
			h = new_point[3];

			function_list->Add(x, v);
			//Печать в таблицу
			dataGridView1->Rows->Add();
			dataGridView1->Rows[i]->Cells[0]->Value = i;
			dataGridView1->Rows[i]->Cells[1]->Value = x;
			dataGridView1->Rows[i]->Cells[2]->Value = v;
			dataGridView1->Rows[i]->Cells[3]->Value = v_2h;
			dataGridView1->Rows[i]->Cells[4]->Value = abs(v - v_2h);
			dataGridView1->Rows[i]->Cells[5]->Value = abs(v - v_2h) * pow(2, p); if (abs(v - v_2h) * pow(2, p) > max_OLP) { max_OLP = abs(v - v_2h) * pow(2, p); max_OLP_x = x; }
			dataGridView1->Rows[i]->Cells[6]->Value = x - last_point[0]; if (x - last_point[0] > max_h) { max_h = x - last_point[0]; max_h_x = x; } if (x - last_point[0] < min_h) { min_h = x - last_point[0]; min_h_x = x; }
			dataGridView1->Rows[i]->Cells[7]->Value = C1;
			dataGridView1->Rows[i]->Cells[8]->Value = C2;
			//dataGridView1->Rows[i]->Cells[9]->Value = 0;
			//dataGridView1->Rows[i]->Cells[10]->Value = 0;
			i++;
			C1_amount += C1;
			C2_amount += C2;
			C1 = 0;
			C2 = 0;
		}
		LineItem^ Curve1 = panel->AddCurve("function_1(x)", function_list, Color::Blue, SymbolType::None);
		//LineItem Curve2 = panel->AddCurve("test_function(x)", test_list, Color::Blue, SymbolType::None);

		textBox12->AppendText(Convert::ToString(i));
		textBox13->AppendText(Convert::ToString(xmax - x));
		textBox14->AppendText(Convert::ToString(max_OLP));
		textBox15->AppendText(Convert::ToString(max_h));
		textBox16->AppendText(Convert::ToString(max_h_x));
		textBox17->AppendText(Convert::ToString(min_h));
		textBox18->AppendText(Convert::ToString(min_h_x));
		textBox19->AppendText(Convert::ToString(max_OLP_x));
		textBox20->AppendText(Convert::ToString(C2_amount));
		textBox21->AppendText(Convert::ToString(C1_amount));
		textBox22->AppendText(Convert::ToString(0));
		textBox23->AppendText(Convert::ToString(0));

		// Устанавливаем интересующий нас интервал по оси X
		panel->XAxis->Scale->Min = xmin_limit;
		panel->XAxis->Scale->Max = xmax_limit;
		/*
				// Устанавливаем интересующий нас интервал по оси Y
				panel->YAxis->Scale->Min = ymin_limit;
				panel->YAxis->Scale->Max = ymax_limit;
		*/
		// Вызываем метод AxisChange (), чтобы обновить данные об осях. 
		// В противном случае на рисунке будет показана только часть графика, 
		// которая умещается в интервалы по осям, установленные по умолчанию
		zedGraphControl1->AxisChange();
		// Обновляем график
		zedGraphControl1->Invalidate();
	}
	private: System::Void zedGraphControl1_Load(System::Object^ sender, System::EventArgs^ e) {};

// Задача 1 без контроля локальной погрешности
private: System::Void button5_Click(System::Object^ sender, System::EventArgs^ e) {

	GraphPane^ panel = zedGraphControl1->GraphPane;
	panel->CurveList->Clear();
	PointPairList^ test_list = gcnew ZedGraph::PointPairList();

	// Интервал, где есть данные
	double xmin = Convert::ToDouble(textBox1->Text);
	double xmax = Convert::ToDouble(textBox2->Text);
	double border = Convert::ToDouble(textBox4->Text);
	double h = Convert::ToDouble(textBox3->Text);
	double u0 = Convert::ToDouble(textBox6->Text);
	std::size_t Max_steps = Convert::ToInt64(textBox11->Text);

	textBox12->Clear();
	textBox13->Clear();
	textBox14->Clear();
	textBox15->Clear();
	textBox16->Clear();
	textBox17->Clear();
	textBox18->Clear();
	textBox19->Clear();
	textBox20->Clear();
	textBox21->Clear();
	textBox22->Clear();
	textBox23->Clear();

	double xmin_limit = xmin - 0.1;
	double xmax_limit = xmax + 0.1;
	// Список точек
	int i = 1;
	dataGridView1->Rows->Clear();
	double v = u0;
	double v_last = u0;
	double v_2h;

	dataGridView1->Rows->Add();
	dataGridView1->Rows[0]->Cells[0]->Value = 0.0;
	dataGridView1->Rows[0]->Cells[1]->Value = xmin;
	dataGridView1->Rows[0]->Cells[2]->Value = floor(v * 10000) / 10000;

	test_list->Add(xmin, v);

	double x = xmin + h;
	double max_OLP = 0.0;
	double max_OLP_x = 0.0;
	double min_h = h;
	double max_h = 0;
	double min_h_x = 0, max_h_x = 0;

	for (; (x <= xmax) && (i < Max_steps); x += h)
	{
		//Добавление на график
		v_last = v;
		v = Runge_Kytta_4(function_1, h, x, v).second;
		test_list->Add(x, v);

		//Печать в таблицу
		v_2h = Runge_Kytta_4(function_1, h / 2.0, x + h / 2.0, Runge_Kytta_4(function_1, h / 2.0, x, v_last).second).second;
		dataGridView1->Rows->Add();
		dataGridView1->Rows[i]->Cells[0]->Value = i;
		dataGridView1->Rows[i]->Cells[1]->Value = x;
		dataGridView1->Rows[i]->Cells[2]->Value = floor(v * 10000) / 10000;
		dataGridView1->Rows[i]->Cells[3]->Value = floor(v_2h * 10000) / 10000;
		dataGridView1->Rows[i]->Cells[4]->Value = abs(v - v_2h);
		dataGridView1->Rows[i]->Cells[5]->Value = S(v, v_2h) * pow(2, p); if (S(v, v_2h) * pow(2, p) > max_OLP) { max_OLP = S(v, v_2h) * pow(2, p); max_OLP_x = x; }
		dataGridView1->Rows[i]->Cells[6]->Value = h; if (h > max_h) { max_h = h; max_h_x = x; } if (h < min_h) { min_h = h; min_h_x = x; }
		dataGridView1->Rows[i]->Cells[7]->Value = 0.0;
		dataGridView1->Rows[i]->Cells[8]->Value = 0.0;
		i++;
	}
	LineItem^ Curve2 = panel->AddCurve("function_1(x)", test_list, Color::Blue, SymbolType::None);

	textBox12->AppendText(Convert::ToString(i));
	textBox13->AppendText(Convert::ToString(xmax - x+h));
	textBox14->AppendText(Convert::ToString(max_OLP));
	textBox15->AppendText(Convert::ToString(max_h));
	textBox16->AppendText(Convert::ToString(max_h_x));
	textBox17->AppendText(Convert::ToString(min_h));
	textBox18->AppendText(Convert::ToString(min_h_x));
	textBox19->AppendText(Convert::ToString(max_OLP_x));
	textBox20->AppendText(Convert::ToString(0));
	textBox21->AppendText(Convert::ToString(0));
	textBox22->AppendText(Convert::ToString(0));
	textBox23->AppendText(Convert::ToString(0));

	// Устанавливаем интересующий нас интервал по оси X
	panel->XAxis->Scale->Min = xmin_limit;
	panel->XAxis->Scale->Max = xmax_limit;
	/*
			// Устанавливаем интересующий нас интервал по оси Y
			panel->YAxis->Scale->Min = ymin_limit;
			panel->YAxis->Scale->Max = ymax_limit;
	*/
	// Вызываем метод AxisChange (), чтобы обновить данные об осях. 
	// В противном случае на рисунке будет показана только часть графика, 
	// которая умещается в интервалы по осям, установленные по умолчанию
	zedGraphControl1->AxisChange();
	// Обновляем график
	zedGraphControl1->Invalidate();

}
// Задача 2 без контроля локальной погрешности
private: System::Void button6_Click(System::Object^ sender, System::EventArgs^ e) {

	GraphPane^ panel = zedGraphControl1->GraphPane;
	panel->CurveList->Clear();
	PointPairList^ function_2_list_1 = gcnew ZedGraph::PointPairList();
	PointPairList^ function_2_list_2 = gcnew ZedGraph::PointPairList();

	// Интервал, где есть данные
	double xmin = Convert::ToDouble(textBox1->Text);
	double xmax = Convert::ToDouble(textBox2->Text);
	double border = Convert::ToDouble(textBox4->Text);
	double h = Convert::ToDouble(textBox3->Text);
	double u0_1 = Convert::ToDouble(textBox5->Text);
	double u0_2 = Convert::ToDouble(textBox8->Text);
	double a = Convert::ToDouble(textBox9->Text);
	double b = Convert::ToDouble(textBox10->Text);
	std::size_t Max_steps = Convert::ToInt64(textBox11->Text);

	textBox12->Clear();
	textBox13->Clear();
	textBox14->Clear();
	textBox15->Clear();
	textBox16->Clear();
	textBox17->Clear();
	textBox18->Clear();
	textBox19->Clear();
	textBox20->Clear();
	textBox21->Clear();
	textBox22->Clear();
	textBox23->Clear();

	double xmin_limit = xmin - 0.1;
	double xmax_limit = xmax + 0.1;
	// Список точек
	int i = 1;
	dataGridView1->Rows->Clear();
	double v_1 = u0_1;
	double v_2 = u0_2;
	double v_last_1 = u0_1;
	double v_last_2 = u0_2;
	double v_2h_1 = 0;
	double v_2h_2 = 0;

	dataGridView1->Rows->Add();
	dataGridView1->Rows[0]->Cells[0]->Value = 0.0;
	dataGridView1->Rows[0]->Cells[1]->Value = xmin;
	dataGridView1->Rows[0]->Cells[2]->Value = floor(v_1 * 10000) / 10000;

	function_2_list_1->Add(xmin, v_1);
	function_2_list_2->Add(xmin, v_2);

	double x = xmin + h;
	double max_OLP = 0.0;
	double max_OLP_x = 0.0;
	double min_h = h;
	double max_h = 0;
	double min_h_x = 0, max_h_x = 0;

	for (; (x <= xmax) && (i < Max_steps); x += h)
	{
		//Добавление на график
		v_last_1 = v_1;
		v_last_2 = v_2;

		v_1 = std::get<1>(Runge_Kytta_4_system(function_2, h, x, v_1, v_2, a, b));
		v_2 = std::get<2>(Runge_Kytta_4_system(function_2, h, x, v_1, v_2, a, b));

		function_2_list_1->Add(x, v_1);
		function_2_list_2->Add(x, v_2);
		//Печать в таблицу
		v_2h_1 = std::get<1>(Runge_Kytta_4_system(function_2, h / 2.0, x + h / 2.0, std::get<1>(Runge_Kytta_4_system(function_2, h / 2.0, x, v_last_1, v_last_2, a, b)), std::get<2>(Runge_Kytta_4_system(function_2, h / 2.0, x, v_last_1, v_last_2, a, b)), a, b));
		dataGridView1->Rows->Add();
		dataGridView1->Rows[i]->Cells[0]->Value = i;
		dataGridView1->Rows[i]->Cells[1]->Value = x;
		dataGridView1->Rows[i]->Cells[2]->Value = floor(v_1 * 10000) / 10000;
		dataGridView1->Rows[i]->Cells[3]->Value = floor(v_2h_1 * 10000) / 10000;
		dataGridView1->Rows[i]->Cells[4]->Value = abs(v_1 - v_2h_1); 
		dataGridView1->Rows[i]->Cells[5]->Value = S(v_1, v_2h_1) * pow(2, p); if (S(v_1, v_2h_1) * pow(2, p) > max_OLP) { max_OLP = S(v_1, v_2h_1) * pow(2, p); max_OLP_x = x; }
		dataGridView1->Rows[i]->Cells[6]->Value = h; if (h > max_h) { max_h = h; max_h_x = x; } if (h < min_h) { min_h = h; min_h_x = x; }
		dataGridView1->Rows[i]->Cells[7]->Value = 0.0;
		dataGridView1->Rows[i]->Cells[8]->Value = 0.0;
		i++;
	}
	LineItem^ Curve1 = panel->AddCurve("v_1(x)", function_2_list_1, Color::Green, SymbolType::None);
	LineItem^ Curve2 = panel->AddCurve("v_2(x)", function_2_list_2, Color::Blue, SymbolType::None);

	textBox12->AppendText(Convert::ToString(i));
	textBox13->AppendText(Convert::ToString(xmax - x+h));
	textBox14->AppendText(Convert::ToString(max_OLP));
	textBox15->AppendText(Convert::ToString(max_h));
	textBox16->AppendText(Convert::ToString(max_h_x));
	textBox17->AppendText(Convert::ToString(min_h));
	textBox18->AppendText(Convert::ToString(min_h_x));
	textBox19->AppendText(Convert::ToString(max_OLP_x));
	textBox20->AppendText(Convert::ToString(0));
	textBox21->AppendText(Convert::ToString(0));
	textBox22->AppendText(Convert::ToString(0));
	textBox23->AppendText(Convert::ToString(0));

	// Устанавливаем интересующий нас интервал по оси X
	panel->XAxis->Scale->Min = xmin_limit;
	panel->XAxis->Scale->Max = xmax_limit;
	/*
			// Устанавливаем интересующий нас интервал по оси Y
			panel->YAxis->Scale->Min = ymin_limit;
			panel->YAxis->Scale->Max = ymax_limit;
	*/
	// Вызываем метод AxisChange (), чтобы обновить данные об осях. 
	// В противном случае на рисунке будет показана только часть графика, 
	// которая умещается в интервалы по осям, установленные по умолчанию
	zedGraphControl1->AxisChange();
	// Обновляем график
	zedGraphControl1->Invalidate();

}
// Задача 2 без контроля локальной погрешности
private: System::Void button7_Click(System::Object^ sender, System::EventArgs^ e) {
	GraphPane^ panel = zedGraphControl1->GraphPane;
	panel->CurveList->Clear();
	PointPairList^ function_1_list = gcnew ZedGraph::PointPairList();
	PointPairList^ function_2_list = gcnew ZedGraph::PointPairList();

	// Интервал, где есть данные
	double xmin = Convert::ToDouble(textBox1->Text);
	double xmax = Convert::ToDouble(textBox2->Text);
	double e_0 = Convert::ToDouble(textBox7->Text);
	double border = Convert::ToDouble(textBox4->Text);
	double u0_1 = Convert::ToDouble(textBox5->Text);
	double u0_2 = Convert::ToDouble(textBox8->Text);
	double a = Convert::ToDouble(textBox9->Text);
	double b = Convert::ToDouble(textBox10->Text);
	std::size_t Max_steps = Convert::ToInt64(textBox11->Text);

	textBox12->Clear();
	textBox13->Clear();
	textBox14->Clear();
	textBox15->Clear();
	textBox16->Clear();
	textBox17->Clear();
	textBox18->Clear();
	textBox19->Clear();
	textBox20->Clear();
	textBox21->Clear();
	textBox22->Clear();
	textBox23->Clear();

	double h = Convert::ToDouble(textBox3->Text);

	double xmin_limit = xmin - 0.1;
	double xmax_limit = xmax + 0.1;

	// Список точек
	int i = 1;
	dataGridView1->Rows->Clear();
	double v_1 = u0_1;
	double v_2 = u0_2;
	double v_last_1 = u0_1;
	double v_last_2 = u0_2;
	double v_2h_1;
	double v_2h_2;
	std::size_t C1 = 0, C2 = 0, last_h = 0;
	std::size_t C1_amount = 0, C2_amount = 0;

	//test_list->Add(xmin, u0);
	std::vector<double> new_point, last_point;
	new_point = { xmin, v_1, v_2, h };

	function_1_list->Add(xmin, v_1);
	function_2_list->Add(xmin, v_2);

	dataGridView1->Rows->Add();
	dataGridView1->Rows[0]->Cells[0]->Value = 0.0;
	dataGridView1->Rows[0]->Cells[1]->Value = xmin;
	dataGridView1->Rows[0]->Cells[2]->Value = floor(u0_1 * 10000) / 10000;

	double x = xmin;
	double max_OLP = 0;
	double max_OLP_x = 0;
	double min_h = h;
	double max_h = 0;
	double min_h_x = 0, max_h_x = 0;

	for (; (x < xmax - border) && (i < Max_steps); ) {
		//Добавление на график
		last_h = h;
		last_point = new_point;

		new_point = RK_4_OLP_for_system(function_2, x, v_1, v_2, h, e_0, a, b);

		x = new_point[0];
		if (x > xmax) {
			new_point = RK_4_OLP_for_system(function_2, last_point[0], v_1, v_2, xmax - last_point[0], e_0, a, b);
			x = new_point[0];
			if (xmax - last_point[0] - h > last_h) {
				C1 = 0; C2 = 0;
			}
		}
		else {
			if (new_point[6] < 0)
				C1 = new_point[6] * -1;
			if (new_point[6] > 0) {
				C2 = 1;
				if (x + h > xmax)
					C2 = 0;
			}
		}

		x = new_point[0];
		v_1 = new_point[1];
		v_2 = new_point[2];
		v_2h_1 = new_point[3];
		v_2h_2 = new_point[4];
		h = new_point[5];

		function_1_list->Add(x, v_1);
		function_2_list->Add(x, v_2);
		//Печать в таблицу
		dataGridView1->Rows->Add();
		dataGridView1->Rows[i]->Cells[0]->Value = i;
		dataGridView1->Rows[i]->Cells[1]->Value = x;
		dataGridView1->Rows[i]->Cells[2]->Value = v_1;
		dataGridView1->Rows[i]->Cells[3]->Value = v_2h_1;
		dataGridView1->Rows[i]->Cells[4]->Value = abs(v_1 - v_2h_1);
		dataGridView1->Rows[i]->Cells[5]->Value = abs(v_1 - v_2h_1) * pow(2, p); if (abs(v_1 - v_2h_1) * pow(2, p) > max_OLP) { max_OLP = abs(v_1 - v_2h_1) * pow(2, p); max_OLP_x = x; }
		dataGridView1->Rows[i]->Cells[6]->Value = x - last_point[0]; if (x - last_point[0] > max_h) { max_h = x - last_point[0]; max_h_x = x; } if (x - last_point[0] < min_h) { min_h = x - last_point[0]; min_h_x = x; }		dataGridView1->Rows[i]->Cells[7]->Value = C1;
		dataGridView1->Rows[i]->Cells[8]->Value = C2;
		//dataGridView1->Rows[i]->Cells[9]->Value = 0;
		//dataGridView1->Rows[i]->Cells[10]->Value = 0;
		C1_amount += C1;
		C2_amount += C2;
		C1 = 0;
		C2 = 0;
		i++;
	}
	LineItem^ Curve1 = panel->AddCurve("v_1(x)", function_1_list, Color::Green, SymbolType::None);
	LineItem^ Curve2 = panel->AddCurve("v_2(x)", function_2_list, Color::Blue, SymbolType::None);

	textBox12->AppendText(Convert::ToString(i));
	textBox13->AppendText(Convert::ToString(xmax - x));
	textBox14->AppendText(Convert::ToString(max_OLP));
	textBox15->AppendText(Convert::ToString(max_h));
	textBox16->AppendText(Convert::ToString(max_h_x));
	textBox17->AppendText(Convert::ToString(min_h));
	textBox18->AppendText(Convert::ToString(min_h_x));
	textBox19->AppendText(Convert::ToString(max_OLP_x));
	textBox20->AppendText(Convert::ToString(C2_amount));
	textBox21->AppendText(Convert::ToString(C1_amount));
	textBox22->AppendText(Convert::ToString(0));
	textBox23->AppendText(Convert::ToString(0));

	// Устанавливаем интересующий нас интервал по оси X
	panel->XAxis->Scale->Min = xmin_limit;
	panel->XAxis->Scale->Max = xmax_limit;
	/*
			// Устанавливаем интересующий нас интервал по оси Y
			panel->YAxis->Scale->Min = ymin_limit;
			panel->YAxis->Scale->Max = ymax_limit;
	*/
	// Вызываем метод AxisChange (), чтобы обновить данные об осях. 
	// В противном случае на рисунке будет показана только часть графика, 
	// которая умещается в интервалы по осям, установленные по умолчанию
	zedGraphControl1->AxisChange();
	// Обновляем график
	zedGraphControl1->Invalidate();

}
};
}
