#define _USE_MATH_DEFINES
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <limits>
#include <sstream>
#include <cctype>
#include <stdexcept>
#include <set>
#include <iomanip>
#include <algorithm>
#include <memory>
#include <unordered_map>
#include <functional>

using namespace std;

// Constants
const double EPSILON = 1e-12;
const double NEWTON_TOLERANCE = 1e-10;
const int MAX_NEWTON_ITERATIONS = 1000;
const double DERIVATIVE_STEP = 1e-8;
const double SOLUTION_TOLERANCE = 1e-8;

// Token types for lexical analysis
enum class TokenType {
    NUMBER, VARIABLE, PLUS, MINUS, MULTIPLY, DIVIDE, POWER,
    LEFT_PAREN, RIGHT_PAREN, FUNCTION, END_OF_INPUT, INVALID
};

struct Token {
    TokenType type;
    string value;
    double number;
    
    Token(TokenType t = TokenType::INVALID, const string& v = "", double n = 0.0)
        : type(t), value(v), number(n) {}
};

// Abstract syntax tree nodes
class ASTNode {
public:
    virtual ~ASTNode() = default;
    virtual double evaluate(double x_value) const = 0;
    virtual string to_string() const = 0;
    virtual unique_ptr<ASTNode> clone() const = 0;
};

class NumberNode : public ASTNode {
private:
    double value;
public:
    NumberNode(double v) : value(v) {}
    double evaluate(double x_value) const override { return value; }
    string to_string() const override { return std::to_string(value); }
    unique_ptr<ASTNode> clone() const override {
        return make_unique<NumberNode>(value);
    }
};

class VariableNode : public ASTNode {
public:
    double evaluate(double x_value) const override { return x_value; }
    string to_string() const override { return "x"; }
    unique_ptr<ASTNode> clone() const override {
        return make_unique<VariableNode>();
    }
};

class BinaryOpNode : public ASTNode {
private:
    unique_ptr<ASTNode> left, right;
    char op;
public:
    BinaryOpNode(unique_ptr<ASTNode> l, char o, unique_ptr<ASTNode> r)
        : left(std::move(l)), op(o), right(std::move(r)) {}
    
    double evaluate(double x_value) const override {
        double l_val = left->evaluate(x_value);
        double r_val = right->evaluate(x_value);
        
        switch (op) {
            case '+': return l_val + r_val;
            case '-': return l_val - r_val;
            case '*': return l_val * r_val;
            case '/': 
                if (abs(r_val) < EPSILON) {
                    return numeric_limits<double>::quiet_NaN();
                }
                return l_val / r_val;
            case '^': return pow(l_val, r_val);
            default: return numeric_limits<double>::quiet_NaN();
        }
    }
    
    string to_string() const override {
        return "(" + left->to_string() + " " + op + " " + right->to_string() + ")";
    }
    
    unique_ptr<ASTNode> clone() const override {
        return make_unique<BinaryOpNode>(left->clone(), op, right->clone());
    }
};

class UnaryOpNode : public ASTNode {
private:
    unique_ptr<ASTNode> operand;
    char op;
public:
    UnaryOpNode(char o, unique_ptr<ASTNode> operand)
        : operand(std::move(operand)), op(o) {}
    
    double evaluate(double x_value) const override {
        double val = operand->evaluate(x_value);
        switch (op) {
            case '-': return -val;
            case '+': return val;
            default: return numeric_limits<double>::quiet_NaN();
        }
    }
    
    string to_string() const override {
        return string(1, op) + operand->to_string();
    }
    
    unique_ptr<ASTNode> clone() const override {
        return make_unique<UnaryOpNode>(op, operand->clone());
    }
};

class FunctionNode : public ASTNode {
private:
    string func_name;
    unique_ptr<ASTNode> argument;
    
    static unordered_map<string, function<double(double)>> functions;
    
public:
    FunctionNode(const string& name, unique_ptr<ASTNode> arg)
        : func_name(name), argument(std::move(arg)) {}
    
    double evaluate(double x_value) const override {
        double arg_val = argument->evaluate(x_value);
        
        auto it = functions.find(func_name);
        if (it != functions.end()) {
            return it->second(arg_val);
        }
        
        return numeric_limits<double>::quiet_NaN();
    }
    
    string to_string() const override {
        return func_name + "(" + argument->to_string() + ")";
    }
    
    unique_ptr<ASTNode> clone() const override {
        return make_unique<FunctionNode>(func_name, argument->clone());
    }
    
    static void initialize_functions() {
        functions["sin"] = [](double x) { return sin(x); };
        functions["cos"] = [](double x) { return cos(x); };
        functions["tan"] = [](double x) { return tan(x); };
        functions["exp"] = [](double x) { return exp(x); };
        functions["log"] = [](double x) { return log(x); };
        functions["ln"] = [](double x) { return log(x); };
        functions["sqrt"] = [](double x) { return sqrt(x); };
        functions["abs"] = [](double x) { return abs(x); };
        functions["asin"] = [](double x) { return asin(x); };
        functions["acos"] = [](double x) { return acos(x); };
        functions["atan"] = [](double x) { return atan(x); };
    }
};

// Static member definition
unordered_map<string, function<double(double)>> FunctionNode::functions;

// Lexical analyzer (tokenizer)
class Lexer {
private:
    string input;
    size_t pos;
    
    void skip_whitespace() {
        while (pos < input.length() && isspace(input[pos])) {
            pos++;
        }
    }
    
    double read_number() {
        size_t start = pos;
        if (pos < input.length() && (input[pos] == '+' || input[pos] == '-')) {
            pos++;
        }
        
        while (pos < input.length() && (isdigit(input[pos]) || input[pos] == '.')) {
            pos++;
        }
        
        if (pos < input.length() && (input[pos] == 'e' || input[pos] == 'E')) {
            pos++;
            if (pos < input.length() && (input[pos] == '+' || input[pos] == '-')) {
                pos++;
            }
            while (pos < input.length() && isdigit(input[pos])) {
                pos++;
            }
        }
        
        string num_str = input.substr(start, pos - start);
        return stod(num_str);
    }
    
    string read_identifier() {
        size_t start = pos;
        while (pos < input.length() && (isalnum(input[pos]) || input[pos] == '_')) {
            pos++;
        }
        return input.substr(start, pos - start);
    }
    
public:
    Lexer(const string& expr) : input(expr), pos(0) {}
    
    Token next_token() {
        skip_whitespace();
        
        if (pos >= input.length()) {
            return Token(TokenType::END_OF_INPUT);
        }
        
        char current = input[pos];
        
        // Numbers
        if (isdigit(current) || current == '.') {
            double num = read_number();
            return Token(TokenType::NUMBER, "", num);
        }
        
        // Identifiers (variables and functions)
        if (isalpha(current) || current == '_') {
            string identifier = read_identifier();
            
            // Check if it's followed by '(' (function)
            skip_whitespace();
            if (pos < input.length() && input[pos] == '(') {
                return Token(TokenType::FUNCTION, identifier);
            }
            
            // Check if it's a variable (x)
            if (identifier == "x") {
                return Token(TokenType::VARIABLE, identifier);
            }
            
            // Check if it's a constant
            if (identifier == "pi" || identifier == "PI") {
                return Token(TokenType::NUMBER, identifier, M_PI);
            }
            if (identifier == "e" || identifier == "E") {
                return Token(TokenType::NUMBER, identifier, M_E);
            }
            
            return Token(TokenType::INVALID, identifier);
        }
        
        // Single character tokens
        pos++;
        switch (current) {
            case '+': return Token(TokenType::PLUS, "+");
            case '-': return Token(TokenType::MINUS, "-");
            case '*': return Token(TokenType::MULTIPLY, "*");
            case '/': return Token(TokenType::DIVIDE, "/");
            case '^': return Token(TokenType::POWER, "^");
            case '(': return Token(TokenType::LEFT_PAREN, "(");
            case ')': return Token(TokenType::RIGHT_PAREN, ")");
            default: return Token(TokenType::INVALID, string(1, current));
        }
    }
};

// Recursive descent parser
class Parser {
private:
    Lexer lexer;
    Token current_token;
    
    void consume(TokenType expected) {
        if (current_token.type != expected) {
            throw runtime_error("Unexpected token: expected " + to_string((int)expected) + 
                              ", got " + to_string((int)current_token.type));
        }
        current_token = lexer.next_token();
    }
    
    // Grammar:
    // expression := term (('+' | '-') term)*
    unique_ptr<ASTNode> parse_expression() {
        auto left = parse_term();
        
        while (current_token.type == TokenType::PLUS || current_token.type == TokenType::MINUS) {
            char op = (current_token.type == TokenType::PLUS) ? '+' : '-';
            consume(current_token.type);
            auto right = parse_term();
            left = make_unique<BinaryOpNode>(std::move(left), op, std::move(right));
        }
        
        return left;
    }
    
    // term := factor (('*' | '/') factor)*
    unique_ptr<ASTNode> parse_term() {
        auto left = parse_factor();
        
        while (current_token.type == TokenType::MULTIPLY || current_token.type == TokenType::DIVIDE) {
            char op = (current_token.type == TokenType::MULTIPLY) ? '*' : '/';
            consume(current_token.type);
            auto right = parse_factor();
            left = make_unique<BinaryOpNode>(std::move(left), op, std::move(right));
        }
        
        return left;
    }
    
    // factor := power | unary
    unique_ptr<ASTNode> parse_factor() {
        return parse_power();
    }
    
    // power := unary ('^' unary)*
    unique_ptr<ASTNode> parse_power() {
        auto left = parse_unary();
        
        if (current_token.type == TokenType::POWER) {
            consume(TokenType::POWER);
            auto right = parse_power(); // Right associative
            left = make_unique<BinaryOpNode>(std::move(left), '^', std::move(right));
        }
        
        return left;
    }
    
    // unary := ('+' | '-') unary | primary
    unique_ptr<ASTNode> parse_unary() {
        if (current_token.type == TokenType::PLUS || current_token.type == TokenType::MINUS) {
            char op = (current_token.type == TokenType::PLUS) ? '+' : '-';
            consume(current_token.type);
            return make_unique<UnaryOpNode>(op, parse_unary());
        }
        
        return parse_primary();
    }
    
    // primary := NUMBER | VARIABLE | FUNCTION '(' expression ')' | '(' expression ')'
    unique_ptr<ASTNode> parse_primary() {
        if (current_token.type == TokenType::NUMBER) {
            double value = current_token.number;
            consume(TokenType::NUMBER);
            return make_unique<NumberNode>(value);
        }
        
        if (current_token.type == TokenType::VARIABLE) {
            consume(TokenType::VARIABLE);
            return make_unique<VariableNode>();
        }
        
        if (current_token.type == TokenType::FUNCTION) {
            string func_name = current_token.value;
            consume(TokenType::FUNCTION);
            consume(TokenType::LEFT_PAREN);
            auto arg = parse_expression();
            consume(TokenType::RIGHT_PAREN);
            return make_unique<FunctionNode>(func_name, std::move(arg));
        }
        
        if (current_token.type == TokenType::LEFT_PAREN) {
            consume(TokenType::LEFT_PAREN);
            auto expr = parse_expression();
            consume(TokenType::RIGHT_PAREN);
            return expr;
        }
        
        throw runtime_error("Unexpected token in primary: " + current_token.value);
    }
    
public:
    Parser(const string& input) : lexer(input) {
        current_token = lexer.next_token();
    }
    
    unique_ptr<ASTNode> parse() {
        auto result = parse_expression();
        if (current_token.type != TokenType::END_OF_INPUT) {
            throw runtime_error("Unexpected input at end: " + current_token.value);
        }
        return result;
    }
};

// Enhanced expression evaluator
class ExpressionEvaluator {
private:
    unique_ptr<ASTNode> ast;
    
public:
    ExpressionEvaluator(const string& expression) {
        try {
            FunctionNode::initialize_functions();
            Parser parser(preprocess_input(expression));
            ast = parser.parse();
        } catch (const exception& e) {
            throw runtime_error("Failed to parse expression '" + expression + "': " + e.what());
        }
    }
    
    double evaluate(double x_value) const {
        if (!ast) {
            return numeric_limits<double>::quiet_NaN();
        }
        return ast->evaluate(x_value);
    }
    
    string to_string() const {
        if (!ast) {
            return "invalid";
        }
        return ast->to_string();
    }
    
    static string preprocess_input(const string& input) {
        string result = input;
        
        // Remove spaces
        result.erase(remove(result.begin(), result.end(), ' '), result.end());
        
        // Add explicit multiplication signs
        for (size_t i = 0; i < result.length() - 1; i++) {
            char current = result[i];
            char next = result[i + 1];
            
            bool need_multiplication = false;
            
            // Number followed by variable or function
            if (isdigit(current) && (next == 'x' || isalpha(next))) {
                need_multiplication = true;
            }
            // ) followed by number, variable, or function
            else if (current == ')' && (isdigit(next) || next == 'x' || isalpha(next))) {
                need_multiplication = true;
            }
            // Variable followed by (
            else if (current == 'x' && next == '(') {
                need_multiplication = true;
            }
            
            if (need_multiplication) {
                result.insert(i + 1, "*");
                i++; // Skip the inserted character
            }
        }
        
        return result;
    }
};

// Numerical derivative
double numerical_derivative(const ExpressionEvaluator& expr, double x) {
    double f1 = expr.evaluate(x + DERIVATIVE_STEP);
    double f2 = expr.evaluate(x - DERIVATIVE_STEP);
    
    if (isnan(f1) || isnan(f2) || isinf(f1) || isinf(f2)) {
        return numeric_limits<double>::quiet_NaN();
    }
    
    return (f1 - f2) / (2 * DERIVATIVE_STEP);
}

// Newton-Raphson method
double newton_method(const ExpressionEvaluator& expr, double initial_guess) {
    double x = initial_guess;
    
    for (int iter = 0; iter < MAX_NEWTON_ITERATIONS; iter++) {
        double fx = expr.evaluate(x);
        
        if (isnan(fx) || isinf(fx)) {
            return numeric_limits<double>::quiet_NaN();
        }
        
        if (abs(fx) < NEWTON_TOLERANCE) {
            return x;
        }
        
        double fpx = numerical_derivative(expr, x);
        if (isnan(fpx) || isinf(fpx) || abs(fpx) < EPSILON) {
            return numeric_limits<double>::quiet_NaN();
        }
        
        double new_x = x - fx / fpx;
        
        if (isnan(new_x) || isinf(new_x)) {
            return numeric_limits<double>::quiet_NaN();
        }
        
        if (abs(new_x - x) < NEWTON_TOLERANCE) {
            return new_x;
        }
        
        x = new_x;
    }
    
    return numeric_limits<double>::quiet_NaN();
}

// Enhanced equation solver
vector<double> solve_equation(const string& expression) {
    set<double> solutions_set;
    vector<double> search_points;
    
    cout << "Processing equation: " << expression << endl;
    
    try {
        ExpressionEvaluator expr(expression);
        cout << "Parsed as: " << expr.to_string() << endl;
        
        // Generate search points
        for (double x = -10.0; x <= 10.0; x += 0.5) {
            search_points.push_back(x);
        }
        
        // Fine search around zero
        for (double x = -2.0; x <= 2.0; x += 0.1) {
            search_points.push_back(x);
        }
        
        // Additional random search points
        for (int i = 0; i < 50; i++) {
            double random_x = -20.0 + (rand() / (double)RAND_MAX) * 40.0;
            search_points.push_back(random_x);
        }
        
        // Apply Newton's method from search points
        cout << "Applying Newton's method from " << search_points.size() << " starting points..." << endl;
        
        for (double start : search_points) {
            double solution = newton_method(expr, start);
            if (!isnan(solution) && !isinf(solution) && abs(solution) < 100) {
                // Verify the solution
                double verification = expr.evaluate(solution);
                if (abs(verification) <= SOLUTION_TOLERANCE) {
                    // Round to avoid floating-point precision issues
                    double rounded = round(solution * 1000000.0) / 1000000.0;
                    solutions_set.insert(rounded);
                }
            }
        }
        
    } catch (const exception& e) {
        cout << "Error: " << e.what() << endl;
        return vector<double>();
    }
    
    vector<double> solutions(solutions_set.begin(), solutions_set.end());
    return solutions;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " \"equation\"" << endl;
        cout << "Examples:" << endl;
        cout << "  " << argv[0] << " \"sin(2*x^4)\"" << endl;
        cout << "  " << argv[0] << " \"x^3 - 2*x + 1\"" << endl;
        cout << "  " << argv[0] << " \"exp(x) - 3\"" << endl;
        cout << "  " << argv[0] << " \"cos(x) + sin(x)\"" << endl;
        return 1;
    }
    
    string input = argv[1];
    
    cout << "Original input: " << input << endl;
    
    vector<double> solutions = solve_equation(input);
    
    if (solutions.empty()) {
        cout << "No real solutions found in the search range." << endl;
    } else {
        cout << "\nFound " << solutions.size() << " solution(s):" << endl;
        for (size_t i = 0; i < solutions.size(); i++) {
            cout << "x" << (i+1) << " = " << fixed << setprecision(6) << solutions[i] << endl;
            
            // Verify each solution
            try {
                ExpressionEvaluator expr(input);
                double check = expr.evaluate(solutions[i]);
                cout << "  Verification: f(" << solutions[i] << ") = " << scientific << setprecision(2) << check << endl;
            } catch (const exception& e) {
                cout << "  Verification failed: " << e.what() << endl;
            }
        }
    }
    
    return 0;
}
