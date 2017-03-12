#include <iostream>

#include <boost/proto/proto.hpp>

using std::cout;
using std::endl;

using boost::proto::_;
using boost::proto::or_;
using boost::proto::plus;
using boost::proto::multiplies;
using boost::proto::terminal;

struct epa : boost::proto::or_<  terminal<_>
                               , plus<epa,epa>
							   , multiplies<epa,epa>
							  >
{ };

int main() {
	terminal<int>::type x{2};
	cout << boost::proto::matches<decltype(x+x), epa>::value << endl;
	cout << boost::proto::matches<decltype(x/!x), epa>::value << endl;
}
