#include <vector>

bool *RLEncode(std::vector<double> s, std::vector<double> &encoded) {
	encoded.clear();

	bool *temp = new bool[s.size()];

	int count = 0;
	int k = 0;

	for (std::vector<double>::size_type i = 0; i < s.size(); i++) {
		if (s[i] == 0 && k > 0) {
			count++;
			if (i == s.size() - 1) {
				encoded.push_back((double)count);
				count = 0;
				temp[k] = true;
				k++;
			}
		}
		else {
			if (count > 0) {
				encoded.push_back((double)count);
				count = 0;
				temp[k] = true;
				k++;
			}
			encoded.push_back(s[i]);
			temp[k] = false;
			k++;
		}
	}
	bool *run = new bool[encoded.size()];

	for (unsigned int i = 0; i < encoded.size(); i++) {
		run[i] = temp[i];
	}

	delete[] temp;

	return run;
}