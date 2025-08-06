#include <iostream>
#include <string>
#include <stack>
#include <vector>
#include <cmath>
#include <limits>
#include <sstream>
#include <cctype>
#include <stdexcept>
#include <tuple>
#include <chrono>
#include <set>
#include <iomanip>
#include <numbers>

using namespace std;

class Monome {
public:
    vector<double> coefficient;
    vector<string> value, ops;

    Monome(vector<double> c, vector<string> v, vector<string> o) : coefficient(c), value(v), ops(o) {}
};

const double rad = 180 / M_PI;
bool deg = 1;

void print(Monome m) {
    cout << "Coeff";
    for(auto i:m.coefficient) cout << 'd' << i;
    cout << '\n';

    cout << "value";
    for(auto i:m.value) cout << 'd' << i;
    cout << '\n';

    cout << "ops";
    for(auto i:m.ops) cout << 'd' << i;
    cout << '\n';
}

Monome compute_output(string str);

string doubleToString(double value) {
    ostringstream oss;
    oss.precision(numeric_limits<double>::digits10 + 1); // Set precision to maximum number of digits for a double
    oss << fixed << value; // Convert double to string
    return oss.str();
}

int find_char(string str, char ch) {
    int output = 0;
    for(auto&i:str) if(i == ch) output++;
    return output;
}

double factorial(double input) {
    if(input < 0) return numeric_limits<double>::infinity();
    return tgamma(input + 1);
}

bool is_number(char ch) {
    return isdigit(ch) || ch == '.' || ch == 'x' || ch == '-';
}

int precedence(char op) {
    if (op == '+' || op == '-') return 1;
    else if (op == '*' || op == '/') return 2;
    else if (op != '(' && op != ')' && op != '[' && op != ']') return 3;
    return 0;
}

double log_base(double base, double number) {
    if(number < 0 || base < 0) return numeric_limits<double>::infinity();
    return log(number) / log(base);
}

bool check_equal_value(string a, string b) {
    vector<string> ch = {"", "(", ")", "[", "]"};
    bool b1 = 0, b2 = 0;
    for(auto& i : ch) {
        if(a == i) b1 = 1;
        if(b == i) b2 = 1;
    }

    return b1 && b2;
}

tuple<double, double> string_to_pow(string a, string b) {
    auto parse_exponent = [](const string& s) -> double {
        if(s == "") return 0.0;
        size_t pos = s.find('^');
        if (pos != string::npos) {
            return stod(s.substr(pos + 1));
        }
        return 1.0; // If no '^' is found, it means it's "x", which is "x^1"
    };

    double exp_a = parse_exponent(a);
    double exp_b = parse_exponent(b);

    return make_tuple(exp_a, exp_b);
}

string remove_parenthesis(string str, bool b) {
    string output;
    int brackets = 0;
    
    for(int i = 0; i < str.size(); i++) {
        if(str[i] == '(' || str[i] == '[') brackets++;
        else if(str[i] == ')' || str[i] == ']') brackets--;

        if(brackets >= 0 || (brackets < 0 && str[i] != ')' && str[i] != ']')) output.push_back(str[i]);
        if(brackets < 0 && (str[i] == ')' || str[i] == ']')) brackets++;
    }

    while(brackets != 0) {
        if(brackets < 0) {output.pop_back(); brackets++;}
        else if(brackets > 0) {output.push_back(')'); brackets--;}
    }
    return output;
}

string cut_all(string str) {
    string output;
    for(auto&i:str) {
        if(i != '(' && i != ')' && i != '[' && i != ']') output.push_back(i);
    }
    return output;
}

string value_doubleToString(tuple<double, double> p, int op) {
    if(get<0>(p) == 0 && get<1>(p) == 0) return "";
    if(op == 0) return "x^" + doubleToString(get<0>(p) - get<1>(p));
    else if(op == 1) return "x^" + doubleToString(get<0>(p) + get<1>(p));
    else if(op == 2) return "x^" + doubleToString(get<0>(p) * get<1>(p));
    return "";
}

// 2x^3 -> x^3

string x_coeff(string str) {
    string output;
    bool b = 0;
    for(auto& i : str) {
        if(i == 'x') b = 1;
        if(b && i != ')' && i != ']') output.push_back(i);
    }
    return output;
}

// string without parenthesis

string cut_par(string str) {
    string output;
    for(auto&i:str) {
        if(i != '(' && i != ')') output.push_back(i);
    }
    return output;
}

// 2x^3 -> 2x

string trim_pow(string str) {
    string out;
    for(auto&i:str) {
        if(i == '^') return out;
        out.push_back(i);
    }
    return out;
}

double coeff(string str) {
    string out;
    bool b = false;
    for(auto&i:str) {
        if(i == '^') b = true;
        else if(b && isdigit(i)) out.push_back(i);
        else if(b && !isdigit(i)) break;
    }

    if(x_coeff(str).size() && out.size() == 0) return 1;
    else if(!x_coeff(str).size() && out.size() == 0) return 0;
    else return stod(out);
}

double coefficient(string str) {
    string s_out;
    for(auto&i:str) {
        if(isdigit(i) || i == '.') s_out.push_back(i);
        else break; 
    }
    return stod(s_out);
}

// 2x^3 -> 3

bool power(string& str1, string& str2) {
    std::string out1, out2;

    bool b = false;
    
    for(auto& i : str1) {
        if(b && ((i < '0' || i > '9') && i != '^' && i != '-' && i != '.' && i != '!' && i != '(' && i != ')')) {
            break;
        } else if(b) {
            out1.push_back(i);
        } else if(i == 'x') {
            b = true; 
            out1.push_back('x');
        }
    }

    b = false;

    for(auto& i : str2) {
        if(b && ((i < '0' || i > '9') && i != '^' && i != '-' && i != '.' && i != '!' && i != '(' && i != ')')) {
            break;
        } else if(b) {
            out2.push_back(i);
        } else if(i == 'x') {
            b = true; 
            out2.push_back('x');
        }
    }
    
    if(find_char(out1, '!') == find_char(out2, '!') && coeff(out1) == coeff(out2)) return 1;
    else return 0;

    if(out1.empty() && out2.empty()) return 1;
    else if(out1.size() == 1 && out2.size() == 1) return 1;
    else if(std::stod(cut_par(out1.substr(2, out1.size()-2))) == std::stod(cut_par(out2.substr(2, out2.size()-2)))) return 1;
    else return 0;
}

// used in exponentiation

string power_coeff(string s1, string s2) {
    string t1, t2;
    bool b = 0;
    for(auto&i : s1) {
        if(i == 'x') b = 1;
        if(b && i != '(' && i != ')' && i != '[' && i != ']') t1.push_back(i);
    }

    b = 0;
    for(auto&i : s2) {
        if(i == 'x') b = 1;
        if(b && i != '(' && i != ')' && i != '[' && i != ']') t2.push_back(i);
    }

    if(t1 == "") t1 = "0";
    if(t2 == "") t2 = s2;
    if(t1 == "x") t1 = "1";
    if(t2 == "x") t2 = "1";
    else if(t1[0] == 'x' && t1[1] == '^') t1 = t1.substr(2, t1.size()-2);
    else if(t2[0] == 'x' && t2[1] == '^') t2 = t2.substr(2, t2.size()-2);
    
    return doubleToString(stod(cut_all(t1))*stod(cut_all(t2)));
}

bool check_product(string a, string b) {
    vector<char> forbidden = {'/'};
    if(a.back() == '!') return 0;
    for(auto&i:a) {
        for(auto&j:forbidden) {
            if (i == j) return 0;
        }
    }
    if(b.back() == '!') return 0;
    for(auto&i:b) {
        for(auto&j:forbidden) {
            if(i == j) return 0;
        }
    }

    return 1;
}

bool check_trig(string a, string b) {
    if(a.size() <= 4 && b.size() <= 4) return 1;
    if(a.size() > 4) {
        for(int i = 0; i < a.size()-3; i++) {
            if(a.substr(i, 3) == "arc" || a.substr(i, 3) == "sin" || a.substr(i, 3) == "cos" || a.substr(i, 3) == "tan" || a.substr(i, 3) == "sec" || a.substr(i, 3) == "csc" || a.substr(i, 3) == "cot" || a.substr(i, 3) == "exp") return 0;
        }
    }
    if(b.size() > 4) {
        for(int i = 0; i < b.size()-3; i++) {
            if(b.substr(i, 3) == "arc" || b.substr(i, 3) == "sin" || b.substr(i, 3) == "cos" || b.substr(i, 3) == "tan" || b.substr(i, 3) == "sec" || b.substr(i, 3) == "csc" || b.substr(i, 3) == "cot" || b.substr(i, 3) == "exp") return 0;
        }
    }
    return 1;
}

string fix_log_brackets(string str) {
    string output;
    int brackets = 0;
    bool comma = 0;

    for(int i = 0; i < str.size(); i++) {
        if(str[i] == '(' || str[i] == '[') brackets++;
        else if(str[i] == ')' || str[i] == ']') brackets--;
    }
    return str;
}

Monome apply_operation(Monome a, Monome b, char op) {
    // SOME PROBLEMS WITH DEFINITION OF INVERSE TRIG FUNCTIONS - OUTPUT IS NOT WHAT DESIRED
    switch(op) {
        case '+': {
            int b1 = a.value.size(), b2 = b.value.size();
            for(int i = 0; i < a.value.size(); i++) {
                for(int j = 0; j < b.value.size(); j++) {
                    if((check_equal_value(x_coeff(a.value[i]), x_coeff(b.value[j])) && power(a.value[i], b.value[j])) && check_trig(a.value[i], b.value[j])) {
                        a.coefficient[i] += b.coefficient[j];
                        b.coefficient.erase(b.coefficient.begin() + j);
                        b.value.erase(b.value.begin() + j);
                        b.ops.erase(b.ops.begin() + j);
                        a.value[i] = "(" + doubleToString(a.coefficient[i]) + cut_par(x_coeff(a.value[i])) + ")";
                        if(b2-- == 0) break;
                    }
                }
            }

            for(auto&i:a.value) i = remove_parenthesis(i,0);
            for(int i = 0; i < b.coefficient.size(); i++) {
                a.coefficient.push_back(b.coefficient[i]);
                a.value.push_back(b.value[i]);
                a.ops.push_back("+");
            }

            if(a.value.size() > 1 && a.value[0][0] != '(' && a.value.back()[a.value.back().size()-1] != ')') {
                a.value[0] = "(" + a.value[0];
                a.value.back() = a.value.back() + ")";
            }
            
            return a;
        }
        case '-': {
            for(int i = 0; i < a.value.size(); i++) {
                for(int j = 0; j < b.value.size(); j++) {
                    if((check_equal_value(x_coeff(a.value[i]), x_coeff(b.value[j])) || power(a.value[i], b.value[j])) && check_trig(a.value[i], b.value[j])) {
                        a.coefficient[i] -= b.coefficient[j];
                        b.coefficient.erase(b.coefficient.begin() + j);
                        b.value.erase(b.value.begin() + j);
                        b.ops.erase(b.ops.begin() + j);
                        a.value[i] = "(" + doubleToString(a.coefficient[i]) + cut_par(x_coeff(a.value[i])) + ")";
                    }
                }
            }

            if(a.ops.back() == "+" && a.coefficient.back() < 0) {
                a.ops.back() = "-";
            }

            for(auto&i:a.value) i = cut_par(i);

            for(int i = 0; i < b.coefficient.size(); i++) {
                a.coefficient.push_back(-b.coefficient[i]);
                a.value.push_back(b.value[i]);
                a.coefficient.back() >= 0? a.ops.push_back("+"): a.ops.push_back("-");
            }
            
            if(a.value.size() > 1 && a.value[0][0] != '(' && a.value.back()[a.value.back().size()-1] != ')') {
                a.value[0] = "(" + a.value[0];
                a.value.back() = a.value.back() + ")";
            }

            return a;
        } 
        case '*': {
            if(check_product(a.value.back(), b.value.back())) {
                Monome m({}, {""}, {""});
                for(int i = 0; i < a.value.size(); i++) {
                    for(int j = 0; j < b.value.size(); j++) {
                        if(find_char(a.value[i], '!') == 0 && find_char(b.value[j], '!') == 0 && check_trig(a.value[i], b.value[j])) {
                            m.coefficient.push_back(a.coefficient[i]*b.coefficient[j]);
                            m.value.push_back(doubleToString(m.coefficient.back()) + value_doubleToString(string_to_pow(x_coeff(a.value[i]), x_coeff(b.value[j])), 1));
                            m.ops.push_back("+");
                        } else if(!check_trig(a.value[i], b.value[j])) {
                            m.coefficient.push_back(1);
                            m.value.push_back(a.value[i] + "*" + b.value[j]);
                            m.ops.push_back("+");
                        } else {
                            m.coefficient.push_back(1);
                            m.value.push_back(doubleToString(a.coefficient[i]) + x_coeff(a.value[i]) + "*" + doubleToString(b.coefficient[j]) + x_coeff(b.value[j]));
                            m.ops.push_back("+");
                        }
                    }
                    
                }
                m.value.erase(m.value.begin());
                for(auto&i:m.value) i = cut_par(i);
                m.value[0] = "(" + m.value[0];
                m.value.back() = m.value.back() + ")";
                m.ops.pop_back();
                return m;
            } else {
                a.ops.push_back("*");
                b.ops.push_back("");
                for(int i = 0; i < b.coefficient.size(); i++) {
                    a.coefficient.push_back(b.coefficient[i]);
                    a.value.push_back(b.value[i]);
                    a.ops.push_back(b.ops[i+1]);
                }
            }

            for(auto&i:a.value) i = cut_par(i);
            a.value[0] = "(" + a.value[0];
            a.value.back() = a.value.back() + ")";

            return a;
        }
        case '/': {
            bool bl = 0;
            for(auto&i:a.value) if(find_char(i, '!') || !check_trig(i, i)) bl = 1;
            for(auto&i:b.value) if(!check_trig(i, i)) bl = 1;
            
            if(b.coefficient.size() > 1 || find_char(b.value.back(), '!') == 1 || bl) {
                a.value[0] = "[" + a.value[0]; a.value.back() += "]"; b.value[0] = "[" + b.value[0]; b.value.back() += "]";
                a.ops.push_back("/");
                b.ops.push_back("");
                for(int i = 0; i < b.value.size(); i++) {
                    a.coefficient.push_back(b.coefficient[i]);
                    a.value.push_back(b.value[i]);
                    a.ops.push_back(b.ops[i+1]);
                }
                //a.value[0] = "[" + a.value[0]; a.value.back() += "]";
                Monome m({0}, {""}, {""});
                for(int i = 0; i < a.value.size(); i++) {
                    m.value.back() += a.value[i];
                    m.value.back() += a.ops[i+1];
                }
                
                //m.value.back() = "[" + remove_parenthesis(m.value.back()) + "]";
                return m;
            }

            for(int i = 0; i < a.value.size(); i++) {
                for(int j = 0; j < b.value.size(); j++) {
                    a.coefficient[i] /= b.coefficient[j];
                    a.value[i] = doubleToString(a.coefficient[i]) + value_doubleToString(string_to_pow(x_coeff(a.value[i]), x_coeff(b.value[j])), 0);
                }
            }
            return a;
        }
        case '^': {
            if(a.coefficient.size() == 1 && b.coefficient.size() == 1 && check_trig(a.value.back(), b.value.back()) && x_coeff(b.value[0]).empty()) {
                if(a.value[0][a.value[0].size()-1] == ')' || (x_coeff(a.value[0]) == "" && x_coeff(b.value[0]) == "")) for(auto&i:a.coefficient) i = pow(i, b.coefficient[0]);
                for(int i = 0; i < a.value.size(); i++) a.value[i] = doubleToString(a.coefficient[i]) + trim_pow(x_coeff(a.value[i])) + "^" + power_coeff(a.value[0], b.value[0]);
                return a; 
            } 

            a.value.back() += ("^[");
            b.ops.push_back("");
            for(int i = 0; i < b.value.size(); i++) {
                a.value.back() += (b.value[i]);
                a.value.back() += (b.ops[i+1]);
            }
            a.value.back() += "]";
            return a;
        }
        case '!': {
            for(int i = 0; i < a.value.size(); i++) {
                if(x_coeff(a.value[0]) == "" && a.coefficient.size() == 1) {
                    a.coefficient[i] = factorial(a.coefficient[i]);
                    a.value[i] = doubleToString(a.coefficient[i]);
                } 
            }

            if(!(x_coeff(a.value[0]) == "" && a.coefficient.size() == 1)) a.value.back().push_back('!');
            return a;
        }
        case 'l': {
            if(a.coefficient.size() == 1 && b.coefficient.size() == 1 && x_coeff(a.value[0]).empty() && x_coeff(b.value[0]).empty()) {
                //cout<<a.coefficient[0]<<' '<<b.coefficient[0];
                a.coefficient[0] = log_base(b.coefficient[0], a.coefficient[0]);
                a.value[0] = doubleToString(a.coefficient[0]);
                return a; 
            }

            Monome m({1}, {""}, {""});
            a.ops.push_back(""); b.ops.push_back("");
            m.value.back() += "log[";

            for(int i = 0; i < a.value.size(); i++) {
                m.value.back() += a.value[i]; m.value.back() += a.ops[i+1];
            }

            m.value.back() += ",";

            for(int i = 0; i < b.value.size(); i++) {
                m.value.back() += b.value[i]; m.value.back() += b.ops[i+1];
            }

            m.value.back() += "]";
            m.value.back() = fix_log_brackets(m.value.back());
            return m;
        }
        case 's': {
            bool bl = 0;
            for(auto&i:a.value) if(!x_coeff(i).empty()) bl = 1; 

            if(bl) {
                a.value[0] = "sin[" + a.value[0];
                a.value.back() = a.value.back() + "]";

                Monome m({1}, {""}, {""});
                a.ops.push_back("");
                for(int i = 0; i < a.value.size(); i++) m.value[0] += a.value[i] + a.ops[i+1];
                return m;
            }

            for(int i = 0; i < a.value.size(); i++) {
                a.coefficient[i] = sin(a.coefficient[i]) * (!deg) + (deg) * sin(a.coefficient[i] / rad);
                a.value[i] = doubleToString(a.coefficient.back());
            }
            return a;
        }
        case 'c': {
            bool bl = 0;
            for(auto&i:a.value) if(!x_coeff(i).empty()) bl = 1;

            if(bl) {
                a.value[0] = "cos[" + a.value[0];
                a.value.back() = a.value.back() + "]";

                Monome m({1}, {""}, {""});
                a.ops.push_back("");
                for(int i = 0; i < a.value.size(); i++) m.value[0] += a.value[i] + a.ops[i+1];
                return m;
            }

            for(int i = 0; i < a.value.size(); i++) {
                a.coefficient[i] = cos(a.coefficient[i]) * (!deg) + (deg) * cos(a.coefficient[i] / rad);
                a.value[i] = doubleToString(a.coefficient.back());
            }
            return a;
        }
        case 't': {
            bool bl = 0;
            for(auto&i:a.value) if(!x_coeff(i).empty()) bl = 1;

            if(bl) {
                a.value[0] = "tan[" + a.value[0];
                a.value.back() = a.value.back() + "]";
                
                Monome m({1}, {""}, {""});
                a.ops.push_back("");
                for(int i = 0; i < a.value.size(); i++) m.value[0] += a.value[i] + a.ops[i+1];
                return m;
            }

            for(int i = 0; i < a.value.size(); i++) {
                a.coefficient[i] = tan(a.coefficient[i]) * (!deg) + (deg) * tan(a.coefficient[i] / rad);
                a.value[i] = doubleToString(a.coefficient.back());
            }
            return a;
        }
        case 'e': {
            if(a.coefficient.size() == 1 && (x_coeff(a.value[0]).empty() || x_coeff(a.value[0])[1] == 'p') && check_trig(a.value[0], a.value[0])) {
                for(int i = 0; i < a.value.size(); i++) {
                    a.coefficient[i] = exp(a.coefficient[i]);
                    a.value[i] = doubleToString(a.coefficient[i]);
                }
                return a;
            } else {
                a.value[0] = "exp[" + a.value[0];
                a.value.back() = a.value.back() + "]";
                return a;
            }
        }
        case 'r': {
            if(a.coefficient.size() == 1 && b.coefficient.size() == 1 && x_coeff(a.value[0]).empty() && x_coeff(b.value[0]).empty()) {
                a.coefficient[0] = pow(a.coefficient[0], 1/b.coefficient[0]);
                a.value[0] = doubleToString(a.coefficient[0]);
                return a; 
            }


            Monome m({1}, {""}, {""});
            a.ops.push_back(""); b.ops.push_back("");
            a.value[0] = "[" + a.value[0];

            for(int i = 0; i < a.value.size(); i++) {
                m.value.back() += a.value[i]; m.value.back() += a.ops[i+1];
            }
            m.value.back() += "]^[1/["; 

            for(int i = 0; i < b.value.size(); i++) {
                m.value.back() += b.value[i]; m.value.back() += b.ops[i+1];
            }

            if(m.value.back().back() == ')') m.value.back().back() = ']';
            m.value.back() += "]]";
            m.value.back() = remove_parenthesis(m.value.back(),0);
            print(m);
            return m;
        }
        case '0': {
            bool bl = 0;
            for(auto&i:a.value) if(!x_coeff(i).empty()) bl = 1;

            if(bl) {
                a.value[0] = "arcsin[" + a.value[0];
                a.value.back() = a.value.back() + "]";
                
                Monome m({1}, {""}, {""});
                a.ops.push_back("");
                for(int i = 0; i < a.value.size(); i++) m.value[0] += a.value[i] + a.ops[i+1];
                return m;
            }

            for(int i = 0; i < a.value.size(); i++) {
                if(isinf(a.coefficient[i])) a.coefficient[i] = numeric_limits<double>::infinity();
                a.coefficient[i] = asin(a.coefficient[i]) * (!deg) + (deg) * rad * asin(a.coefficient[i]);
                a.value[i] = doubleToString(a.coefficient[i]);
            }

            return a;
        }
        case '1': {
            bool bl = 0;
            for(auto&i:a.value) if(!x_coeff(i).empty()) bl = 1;

            if(bl) {
                a.value[0] = "arccos[" + a.value[0];
                a.value.back() = a.value.back() + "]";
                
                Monome m({1}, {""}, {""});
                a.ops.push_back("");
                for(int i = 0; i < a.value.size(); i++) m.value[0] += a.value[i] + a.ops[i+1];
                return m;
            }

            for(int i = 0; i < a.value.size(); i++) {
                a.coefficient[i] = acos(a.coefficient[i]) * (!deg) + (deg) * rad * acos(a.coefficient[i]);
                a.value[i] = doubleToString(a.coefficient[i]);
            }
            return a;
        }
        case '2': {
            bool bl = 0;
            for(auto&i:a.value) if(!x_coeff(i).empty()) bl = 1;

            if(bl) {
                a.value[0] = "arctan[" + a.value[0];
                a.value.back() = a.value.back() + "]";
                
                Monome m({1}, {""}, {""});
                a.ops.push_back("");
                for(int i = 0; i < a.value.size(); i++) m.value[0] += a.value[i] + a.ops[i+1];
                return m;
            }

            for(int i = 0; i < a.value.size(); i++) {
                a.coefficient[i] = atan(a.coefficient[i]) * (!deg) + (deg) * rad * atan(a.coefficient[i]);
                a.value[i] = doubleToString(a.coefficient[i]);
            }
            return a;
        }
        case '3': {
            bool bl = 0;
            for(auto&i:a.value) if(!x_coeff(i).empty()) bl = 1;

            if(bl) {
                a.value[0] = "arccsc[" + a.value[0];
                a.value.back() = a.value.back() + "]";
                
                Monome m({1}, {""}, {""});
                a.ops.push_back("");
                for(int i = 0; i < a.value.size(); i++) m.value[0] += a.value[i] + a.ops[i+1];
                return m;
            }

            for(int i = 0; i < a.value.size(); i++) {
                a.coefficient[i] = asin(1/a.coefficient[i]) * (!deg) + (deg) * rad * asin(1/a.coefficient[i]);
                a.value[i] = doubleToString(a.coefficient[i]);
            }
            return a;
        }
        case '4': {
            bool bl = 0;
            for(auto&i:a.value) if(!x_coeff(i).empty()) bl = 1;

            if(bl) {
                a.value[0] = "arccot[" + a.value[0];
                a.value.back() = a.value.back() + "]";
                
                Monome m({1}, {""}, {""});
                a.ops.push_back("");
                for(int i = 0; i < a.value.size(); i++) m.value[0] += a.value[i] + a.ops[i+1];
                return m;
            }

            for(int i = 0; i < a.value.size(); i++) {
                a.coefficient[i] = acos(1/a.coefficient[i]) * (!deg) + (deg) * rad * acos(1/a.coefficient[i]);
                a.value[i] = doubleToString(a.coefficient[i]);
            }
            return a;
        }
        case '5': {
            bool bl = 0;
            for(auto&i:a.value) if(!x_coeff(i).empty()) bl = 1;

            if(bl) {
                a.value[0] = "arcsec[" + a.value[0];
                a.value.back() = a.value.back() + "]";
                
                Monome m({1}, {""}, {""});
                a.ops.push_back("");
                for(int i = 0; i < a.value.size(); i++) m.value[0] += a.value[i] + a.ops[i+1];
                return m;
            }

            for(int i = 0; i < a.value.size(); i++) {
                a.coefficient[i] = atan(1/a.coefficient[i]) * (!deg) + (deg) * rad * atan(1/a.coefficient[i]);
                a.value[i] = doubleToString(a.coefficient[i]);
            }
            return a;
        }
        case 'w': {
            bool bl = 0;
            for(auto&i:a.value) if(!x_coeff(i).empty()) bl = 1;

            if(bl) {
                a.value[0] = "sec[" + a.value[0];
                a.value.back() = a.value.back() + "]";
                
                Monome m({1}, {""}, {""});
                a.ops.push_back("");
                for(int i = 0; i < a.value.size(); i++) m.value[0] += a.value[i] + a.ops[i+1];
                return m;
            }

            for(int i = 0; i < a.value.size(); i++) {
                a.coefficient[i] = 1/cos(a.coefficient[i]) * (!deg) + (deg) * 1/cos(a.coefficient[i] / rad);
                a.value[i] = doubleToString(a.coefficient.back());
            }
            return a;
        }
        case 'x': {
            bool bl = 0;
            for(auto&i:a.value) if(!x_coeff(i).empty()) bl = 1;

            if(bl) {
                a.value[0] = "csc[" + a.value[0];
                a.value.back() = a.value.back() + "]";
                
                Monome m({1}, {""}, {""});
                a.ops.push_back("");
                for(int i = 0; i < a.value.size(); i++) m.value[0] += a.value[i] + a.ops[i+1];
                return m;
            }

            for(int i = 0; i < a.value.size(); i++) {
                a.coefficient[i] = 1/sin(a.coefficient[i]) * (!deg) + (deg) * 1/sin(a.coefficient[i] / rad);
                a.value[i] = doubleToString(a.coefficient.back());
            }
            return a;
        }
        case 'y': {
            bool bl = 0;
            for(auto&i:a.value) if(!x_coeff(i).empty()) bl = 1;

            if(bl) {
                a.value[0] = "cot[" + a.value[0];
                a.value.back() = a.value.back() + "]";
                
                Monome m({1}, {""}, {""});
                a.ops.push_back("");
                for(int i = 0; i < a.value.size(); i++) m.value[0] += a.value[i] + a.ops[i+1];
                return m;
            }

            for(int i = 0; i < a.value.size(); i++) {
                a.coefficient[i] = 1/tan(a.coefficient[i]) * (!deg) + (deg) * 1/tan(a.coefficient[i] / rad);
                a.value[i] = doubleToString(a.coefficient.back());
            }
            return a;
        }
        default: {
            return Monome({}, {""}, {""}); // default case, though it should never reach here
        }
    }
}

string process_string(string input) {
    string s1, s2;
    bool b = 0;
    int brackets = 0, total = 1;

    for(int i = 0; i < input.size(); i++) {
        if(input[i] == '(' || input[i] == '[') brackets++;
        else if(input[i] == ')' || input[i] == ']') brackets--;
        else if(input[i] == 'r' || input[i] == 'l') total++;
        else if(input[i] == ',') {
            total--;
            if(total == 0) {
                b = 1;
                continue;
            }
        }
        if(!b) s1.push_back(input[i]);
        if(b) s2.push_back(input[i]);
    }

    Monome out1 = compute_output(remove_parenthesis(s1.substr(1, s1.size()-1), 0));
    Monome out2 = compute_output(remove_parenthesis(s2.substr(0, s2.size()-1), 0));

    string output;
    brackets = 0;

    for(int i = 0; i < out1.coefficient.size(); i++) {
       output += out1.ops[i]; output += out1.value[i];
    }

    for(auto&i:output) if(i == '(') {brackets++;} else if(i == ')') {brackets--;}
    while(brackets < 1 && (output.back() == ')' || output.back() == ']')) {output.pop_back(); brackets++;}

    output += ",";

    for(int i = 0; i < out2.coefficient.size(); i++) {
       output += out2.ops[i]; output += out2.value[i];
    }

    return remove_parenthesis(output, 0);
} 

string process_trig(string input) {
    string s1;
    int brackets = 1;
    for(int i = 0; i < input.size(); i++) {
        if(brackets == 1 && input[i] == ')' && input[i] != ']') break;
        if(input[i] == '(' || input[i] == '[') brackets++;
        if(input[i] == ')' || input[i] == ']') brackets--;

        s1.push_back(input[i]);
    }
    for(auto&i:s1) if(i == 'x') return remove_parenthesis(s1.substr(1, s1.size()-2),0);
    return remove_parenthesis(compute_output(s1.substr(1, s1.size()-2)).value.back(),0);
}

void process_operator(vector<Monome>& values, vector<char>& ops) {
    vector<char> two_terms = {'+', '-', '*', '/', '^', 'r', 'l', '(', ')', '[', ']'};
    bool present = 0;
    for(auto& i:two_terms) if(ops[ops.size()-1] == i) present = 1;
    if(present) {
        Monome b = values[values.size()-1]; values.pop_back();
        Monome a = values[values.size()-1]; values.pop_back();
        char op = ops[ops.size()-1]; ops.pop_back();
    
        values.push_back({apply_operation(a, b, op)});
    } else {
        Monome a = values[values.size()-1]; values.pop_back();
        char op = ops[ops.size()-1]; ops.pop_back();
        values.push_back({apply_operation(a, Monome({}, {""}, {""}), op)});
    }
}

string fix_string(string str) {
    // general patches
    string output;
    int brackets = 0;

    if(str[0] == '-') {output = "0+-";}
    else if(str[0] == 'x') output = "1x";
    else output.push_back(str[0]);

    for(int i = 1; i < str.size(); i++) {
        if(!is_number(str[i]) && brackets) {output.push_back(')'); brackets--;}
        if(str[i] == '-' && !is_number(str[i-1]) && str[i-1] != ')') {output += "(0+-"; brackets++;}
        else if(str[i] == '-' && ((is_number(str[i-1]) && str[i-1] != '^') || str[i-1] == ')')) output += "+-";
        else if(str[i] == 'x' && !is_number(str[i-1]) && str[i-1] != 'e') output += "1x";
        else output.push_back(str[i]);
    }
    
    while(brackets) {output += ")"; brackets--;}
    return output;
}

Monome compute_output(string str) {
    vector<Monome> values;
    vector<char> ops;
    string pb_str;
    int br = 0; // brackets

    // lets the user choose the measure unity for angles
    if(str == "mode=deg") {
        deg = 1;
        cout << "Measure unity for angles: degrees" << '\n';
        return Monome({2}, {"deg"}, {""});
    } else if(str == "mode=rad") {
        deg = 0;
        cout << "Measure unity for angles: radians" << '\n';
        return Monome({2}, {"rad"}, {""});
    }

    str = fix_string(str); // string adjusted for bugs

    for (int i = 0; i < str.size(); i++) {

        if (str[i] == ' ') continue;

        else if (is_number(str[i])) pb_str.push_back(str[i]);

        else if (str[i] == 'l' && str.substr(i, 3) == "log" || str[i] == 'r' && str.substr(i, 4) == "root") {
            str[i] == 'r' ? ops.push_back('r') : ops.push_back('l');
            i += 2 + int(str[i] == 'r');  
            string replace;
            int brackets = 0, j = i+1;
            for(j = i+1; j < str.size(); j++) {
                if(str[j] == ')' || str[j] == ']') brackets--;
                if(str[j] == '(' || str[j] == '[') brackets++;

                replace.push_back(str[j]);
                if(brackets == 0) break;
            }
            string repl_str = str.substr(0, i+2) + process_string(replace) + str.substr(j, str.size()-j+1);
            str.clear();
            str = repl_str;

        } else if(str.substr(i, 3) == "sin" || str.substr(i, 3) == "cos" || str.substr(i, 3) == "tan" || str.substr(i, 3) == "sec" || str.substr(i, 3) == "csc" || str.substr(i, 3) == "cot" || str.substr(i, 3) == "exp") {
            if(str.substr(i, 3) == "sin") ops.push_back('s');
            else if(str.substr(i, 3) == "cos") ops.push_back('c');
            else if(str.substr(i, 3) == "tan") ops.push_back('t');
            else if(str.substr(i, 3) == "sec") ops.push_back('w');
            else if(str.substr(i, 3) == "csc") ops.push_back('x');
            else if(str.substr(i, 3) == "cot") ops.push_back('y');
            else if(str.substr(i, 3) == "exp") ops.push_back('e');
            i += 2;

            string replace;
            int brackets = 0, j = i+1;
            for(j = i+1; j < str.size(); j++) {
                if(str[j] == ')' || str[j] == ']') brackets--;
                if(str[j] == '(' || str[j] == '[') brackets++;

                replace.push_back(str[j]);
                if(brackets == 0) break;
            }
            
            string repl_str = str.substr(0, i+2) + process_trig(replace) + str.substr(j, str.size()-j+1);
            str.clear();
            str = repl_str;

        } else if(str.substr(i, 6) == "arcsin" || str.substr(i, 6) == "arccos" || str.substr(i, 6) == "arctan" || str.substr(i, 6) == "arcsec" || str.substr(i, 6) == "arccsc" || str.substr(i, 6) == "arccot") {
            if(str.substr(i, 6) == "arcsin") ops.push_back('0');
            else if(str.substr(i, 6) == "arccos") ops.push_back('1');
            else if(str.substr(i, 6) == "arctan") ops.push_back('2');
            else if(str.substr(i, 6) == "arccsc") ops.push_back('3');
            else if(str.substr(i, 6) == "arccot") ops.push_back('4');
            else if(str.substr(i, 6) == "arcsec") ops.push_back('5');
            i += 5;
            
            string replace;
            int brackets = 0, j = i+1;
            for(j = i+1; j < str.size(); j++) {
                if(str[j] == ')' || str[j] == ']') brackets--;
                if(str[j] == '(' || str[j] == '[') brackets++;

                replace.push_back(str[j]);
                if(brackets == 0) break;
            }
            string repl_str = str.substr(0, i+2) + process_trig(replace) + str.substr(j, str.size()-j+1);
            str.clear();
            str = repl_str;
        } else {
            if (!pb_str.empty() && pb_str[pb_str.size()-1] != '(' && pb_str[pb_str.size()-1] != ')' && pb_str[pb_str.size()-1] != '[' && pb_str[pb_str.size()-1] != ']') { // AND THERE EXISTS A CHAR IN PB_STR THAT IS DIFFERENT FROM ( OR )
                values.push_back(Monome({stod(cut_all(pb_str))}, {pb_str}, {""}));
                pb_str.clear();
            }
            else if((pb_str == "(" || pb_str == "[") && values.size() > 0) {
                pb_str == "(" ? values.back().value.back() = "(" + values.back().value.back() : values.back().value.back() = "[" + values.back().value.back();
                pb_str.clear();
            }
            else if((pb_str == ")" || pb_str == "]") && values.size() > 0) {
                pb_str == ")" ? values.back().value.back() = values.back().value.back() + ")" : values.back().value.back() = values.back().value.back() + "]";
                pb_str.clear();
            }
            if (str[i] == '(' || str[i] == '[') {
                br++;
                str[i] == '(' ? pb_str.push_back('(') : pb_str.push_back('[');
                ops.push_back(str[i]);
            }
            else if (str[i] == ')' || str[i] == ']') {
                br--;
                str[i] == ')' ? pb_str.push_back(')') : pb_str.push_back(']');
                if(br >= 0) while (!ops.empty() && ops.back() != '(' && ops.back() != '[') process_operator(values, ops);
                ops.pop_back();

            } else if (str[i] == ',') {
                while(ops.back() != 'r' && ops.back() != 'l' && ops.back() != '(' && ops.back() != '[') process_operator(values, ops);
            } else {
                while (!ops.empty() && precedence(ops.back()) >= precedence(str[i])) process_operator(values, ops);
                ops.push_back(str[i]);
            }
        }
    }   
    
    if (!pb_str.empty() && pb_str != "(" && pb_str != ")" && pb_str != "[" && pb_str != "]") {
        values.push_back(Monome({stod(cut_all(pb_str))}, {pb_str}, {""}));
    }
    else if(pb_str == "(" || pb_str == "[") {
        pb_str == "(" ? values.back().value.back() = "(" + values.back().value.back() : values.back().value.back() = "[" + values.back().value.back();
        pb_str.clear();
    }
    else if(pb_str == ")" || pb_str == "]") {
        pb_str == ")" ? values.back().value.back() = values.back().value.back() += ")" : values.back().value.back() = values.back().value.back() + "]";
        pb_str.clear();
    }

    while (!ops.empty()) {
        process_operator(values, ops);
    }
    
    while(values.back().coefficient.empty() && values.back().value.empty()) values.pop_back();
    return values.back();
}

string clean_string(string str) {
    string output;
    for(auto i:str) {
        if(is_number(i)) output += i;
    }
    return output;
}

double f(string s, double x) {
    string str; double result;
    if(isinf(x) || isnan(x)) return numeric_limits<double>::infinity();

    if(s[0] != 'x') str.push_back(s[0]);
    else str += "*" + doubleToString(x);

    for(int i = 1; i < s.size(); i++) {
        if(s[i] != 'x' || s[i-1] == 'e') str.push_back(s[i]);
        else str += "*" + doubleToString(x);
    }
    
    str = compute_output(fix_string(str)).value.back();
    if(str.substr(0, 3) == "nan") return numeric_limits<double>::infinity();
    str = clean_string(str);
    
    try {
        result = stod(str);
    } catch (const std::invalid_argument& e) {
        return numeric_limits<double>::infinity();
    } catch (const std::out_of_range& e) {
        return numeric_limits<double>::infinity();
    }

    return result;
}

double f_prime(string s, double x) {
    return (f(s, x + 10e-7) - f(s, x)) / 10e-7;
}

double newton(string s, double x) {
    int iter = -1; double back = x;
    while(abs(f(s, x)) > 10e-15 && iter++ < 2150) {
        x -= 0.5 * f(s, x) / f_prime(s, x);
        if(isinf(x) || isnan(x)) return numeric_limits<double>::infinity();
        else back = x;
    }
    return x;
}

vector<double> solve_eq(string str) {
    set<double> sols;
    vector<double> solutions;

    int start = -50, end = 50;
    double tolerance = 10e-15;

    while (start <= end) {
        double x1 = newton(str, start);
        double x2 = newton(str, end);


        if (std::isinf(x1) && std::isinf(x2)) {
            start++;
            end--;
            continue;
        } else if (std::isinf(x1)) {
            if (std::abs(f(str, x2)) <= tolerance) {
                sols.insert(round(x2 * 1000.0) / 1000.0);
            }
            start++;
            end--;
        } else if (std::isinf(x2)) {
            if (std::abs(f(str, x1)) <= tolerance) {
                sols.insert(round(x1 * 1000.0) / 1000.0);
            }
            start++;
            end--;
        } else {
            if (std::abs(f(str, x1)) <= tolerance) {
                sols.insert(round(x1 * 1000.0) / 1000.0);
            }
            if (std::abs(f(str, x2)) <= tolerance) {
                sols.insert(round(x2 * 1000.0) / 1000.0);
            }

            if (x1 == x2) {
                break;
            } else {
                double mid = (start + end) / 2;
                double x_mid = newton(str, mid);
                
                if (std::abs(f(str, x_mid)) <= tolerance) {
                    sols.insert(round(x_mid * 1000.0) / 1000.0);
                }

                // Adjust the search range
                if (x_mid > start) {
                    end = mid - 1;
                } else {
                    start = mid + 1;
                }
            }
        }
    }
    cout<<"\n\n";
    for(auto i:sols) solutions.push_back(i);
    return solutions;
}

void calculator(string input) {
        string str = input;


        int brackets = 0;
        for(auto& i : str) {
            if(i == '(') brackets++;
            if(i == ')') brackets--;
        }

        if(brackets != 0) {
            throw invalid_argument("Brackets do not match");
        }

        Monome output = compute_output(str);
        
        string result;
        for(int i = 0; i < output.coefficient.size(); i++) {
            result += output.ops[i]; result += output.value[i];
        }
        
        result = remove_parenthesis(result,0);
        if(result != "deg" && result != "rad") {
            vector<double> solutions = solve_eq(result);
            for(auto i:solutions) cout << i << '\n';
        }
    }


int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Error: no equation received.\n";
        return 1;
    }

    std::string input = argv[1];
    calculator(input);
    return 0;
}