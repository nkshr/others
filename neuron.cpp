#include <iostream>
#include <vector>

using namespace std;

struct data
{
	double input[2];
	double true_data;
};

class neuron
{
private:
	double w[2];
	double out[4];
	double theta = 0.5, epsilon = 0.5;
	data data_list[4];
	int counter = 0;

	void calc_output()
	{
		for (int i = 0; i < 4; i++)
		{
			out[i] = 0;
			out[i] += data_list[i].input[0] * w[0] ;
			out[i] += data_list[i].input[1] * w[1] ;
			out[i] -= theta;
			if (out[i] <= 0)
				out[i] = 0;
			else
				out[i] = 1;
		}
	}

	void show_result()
	{
		cout << "///////step" << ++counter
			<< "////////" << endl;
		cout << "input		output" << endl;
		for (int i = 0; i < 4; i++)
		{
			cout << data_list[i].input[0] << " " << data_list[i].input[1];
			cout << "		" << out[i] << endl;
		}
		cout << endl;
	}

	bool check_result()
	{
		for (int i = 0; i < 4; i++)
		{
			if (data_list[i].true_data != out[i])
				return false;
		}
		return true;
	}

	void learn()
	{
		double delta_w;
		for (int i = 0; i < 4; i++)
		{
			delta_w = (-1) * epsilon *
				(out[i] - data_list[i].true_data) * data_list[i].input[0];
			w[0] += delta_w;
			delta_w = (-1) * epsilon *
				(out[i] - data_list[i].true_data) * data_list[i].input[1];
			w[1] += delta_w;
		}
	}


public:
	neuron()
	{
		for (int i = 0; i < 4; i++)
		{
			w[i] = 0;
		}
	};
	
	void add_data(data arg)
	{
		data_list[counter] = arg;
		counter++;
		if (counter == 4)
			counter = 0;
	};

	void run()
	{
		for (;;)
		{
			calc_output();
			show_result();
			if (check_result())
				break;
			learn();
		}
		//cout << w[0] << "	" << w[1] << endl;
	};
};
int main(int argc, char **argv)
{
	/* AND */
	cout << "************AND neuron************" << endl;
	data d1, d2, d3, d4;
	d1.input[0] = 0;
	d1.input[1] = 0;
	d1.true_data = 0;
	
	d2.input[0] = 1;
	d2.input[1] = 0;
	d2.true_data = 0;

	d3.input[0] = 0;
	d3.input[1] = 1;
	d3.true_data = 0;

	d4.input[0] = 1;
	d4.input[1] = 1;
	d4.true_data = 1;
	
	neuron and_neu;
	and_neu.add_data(d1);
	and_neu.add_data(d2);
	and_neu.add_data(d3);
	and_neu.add_data(d4);
	and_neu.run();

	/* OR */
	cout << "**************OR neuron*************" << endl;
	data d5, d6, d7, d8;
	d5.input[0] = 0;
	d5.input[1] = 0;
	d5.true_data = 0;

	d6.input[0] = 1;
	d6.input[1] = 0;
	d6.true_data = 1;

	d7.input[0] = 0;
	d7.input[1] = 1;
	d7.true_data = 1;

	d8.input[0] = 1;
	d8.input[1] = 1;
	d8.true_data = 1;

	neuron or_neu;
	or_neu.add_data(d5);
	or_neu.add_data(d6);
	or_neu.add_data(d7);
	or_neu.add_data(d8);
	or_neu.run();

	/* XOR */
	/*cout << "**************OR neuron*************" << endl;
	data d9, d10, d11, d12;
	d9.input[0] = 0;
	d9.input[1] = 0;
	d9.true_data = 0;

	d10.input[0] = 1;
	d10.input[1] = 0;
	d10.true_data = 1;

	d11.input[0] = 0;
	d11.input[1] = 1;
	d11.true_data = 1;

	d12.input[0] = 1;
	d12.input[1] = 1;
	d12.true_data = 0;

	neuron xor_neu;
	xor_neu.add_data(d9);
	xor_neu.add_data(d10);
	xor_neu.add_data(d11);
	xor_neu.add_data(d12);
	xor_neu.run();*/
	return 0;
}