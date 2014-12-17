#ifndef inputReader_H
#define inputReader_H

#include <cstdlib>
#include <vector>

std::vector<int> inputRead (void)
{
	char n[100];
	std::cin.getline(n,100,'\n');
	std::vector<int> conditions;

	int condition =0;

	for(int i=0; i<100; i++)
	{
		if(n[i] >='0' && n[i] <= '9')
		{
			condition = condition *10;
			condition += n[i] -'0';
		}
		if(n[i]== ',')
		{
			conditions.push_back(condition);		
			std::cout << '\n' << "condition selected:  " << condition;
			condition = 0;
		}
	}
	return conditions;
}
#endif