#include <iostream>

using std::cout;
using std::endl;

#include <boost/proto/proto.hpp>

using boost::proto::_;
using boost::proto::or_;
using boost::proto::when;
using boost::proto::_value;
using boost::proto::_left;
using boost::proto::_right;
using boost::proto::otherwise;
using boost::proto::_default;
using boost::proto::terminal;
using boost::proto::plus;

struct custom : boost::proto::callable {
	using result_type = double;

	template <typename T, typename U>
	result_type operator()(const T& t, const U& u) const {
		return (t+u)*100;
	}
};

struct eval : boost::proto::or_< when< terminal<_>
                                     , _value>
                               , when< plus<eval,eval>
							         , custom(eval(_left), eval(_right)) >
                               , otherwise< _default<eval> >
							   >
{ };

int main() {
	terminal<int>::type x{2};
	eval e;
	cout << e((x+3)+(x*3)) << endl;
}	
