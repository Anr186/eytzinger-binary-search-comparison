#include <iostream>
#include <chrono>
#include <ctime>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <bitset>
#include <string>
#include <memory>
#include <fstream>
#include <vector>
using namespace std;

template<typename T, typename I> 
class sorted_array {
public:
    sorted_array(const T arr[], I size) : a(arr), n(size) {}

    I branchy_search(T x) const {
        I lo = 0;
        I hi = n;
        while (lo < hi) {
            I m = (lo + hi) / 2;
            if (x < a[m]) {
                hi = m;
            } else if (x > a[m]) {
                lo = m + 1;
            } else {
                return m;
            }
        }
        return hi;
    }
    
    I _branchfree_search(T x) const {
        const T *base = a; 
        I n = this->n; 
        while (n > 1) {
            const I half = n / 2;
            base = (base[half] < x) ? &base[half] : base; 
            n -= half;
        }
        return (*base < x) + base - a;
    }
    I _branchfreepref_search(T x) const {
        const T *base = a;
        I n = this->n;
        
        while (n > 1) {
            const I half = n / 2;
            
            
            __builtin_prefetch(&base[half/2], 0, 0);    
            __builtin_prefetch(&base[half + half/2], 0, 0); 
            
            base = (base[half] < x) ? &base[half] : base;
            n -= half;
        }
        
        return (*base < x) + base - a;
    }
    I eytz_branchy_search(T x) const {
        I i = 0;
        while (i < n) {
            if (x < a[i]) {
                i = 2*i + 1;
            }
            else if (x > a[i]) {
                i = 2*i + 2;
            } 
            else {
                return i;
            }
        }
        I j = (i+1) >> __builtin_ffs(~(i+1));
        return (j == 0) ? n : j-1;
    }
    I eytz_branchfree_search(T x) const {
        I i = 0;
        while (i < n) {
            i = (x <= a[i]) ? (2*i + 1) : (2*i + 2);
        }
        I j = (i+1) >> __builtin_ffs(~(i+1));
        return (j == 0) ? n : j-1;
    }
	void print_array() const {
        cout << "Array contents: [";
        for (I i = 0; i < min(n, 50); ++i) { 
            cout << a[i];
            if (i != min(n, 50) - 1) cout << ", ";
        }
        if (n > 50) cout << ", ...";
        cout << "]" << endl;
    }
private:
    const T* a;
    I n;
};

struct TreeNode {
    int data;
    unique_ptr<TreeNode> left;
    unique_ptr<TreeNode> right;
    
    TreeNode(int val) : data(val), left(nullptr), right(nullptr) {}
};

unique_ptr<TreeNode> BuildCompleteBST(vector<int>& arr, int start, int end) {
    if (start > end) return nullptr;

    int mid = start + (end - start) / 2;
    auto node = make_unique<TreeNode>(arr[mid]);
    node->left = BuildCompleteBST(arr, start, mid - 1);
    node->right = BuildCompleteBST(arr, mid + 1, end);
    return node;
}

vector<int> LevelOrderToArray(const unique_ptr<TreeNode>& root) {
    vector<int> result;
    if (!root) return result;

    vector<const TreeNode*> queue;
    queue.push_back(root.get());

    size_t front = 0;
    while (front < queue.size()) {
        const TreeNode* current = queue[front++];
        result.push_back(current->data);

        if (current->left) queue.push_back(current->left.get());
        if (current->right) queue.push_back(current->right.get());
    }

    return result;
}

int pow2(int n) { 
    return 1 << n;
}

void fillmas(int *a, int n) {
    for(int i = 0; i < n; i++) {
        a[i] = rand() % n;
    }
    sort(a, a + n);
}

void fillvec(vector<int>& a, int n) {
    a.resize(n);
    for (int i = 0; i < n; i++) {
        a[i] = rand() % n;
    }
    sort(a.begin(), a.end());
}

void create_bin_files(double avrg_time[], double first_time[], double last_time[], int size, string type = "") {
    ofstream outFile1("AVG_arr_" + type + ".bin", ios::binary);  
    outFile1.write(reinterpret_cast<char*>(avrg_time), size * sizeof(double));
    outFile1.close();

    ofstream outFile2("FRST_arr_" + type + ".bin", ios::binary);  
    outFile2.write(reinterpret_cast<char*>(first_time), size * sizeof(double));
    outFile2.close();

    ofstream outFile3("LAST_arr_" + type + ".bin", ios::binary);  
    outFile3.write(reinterpret_cast<char*>(last_time), size * sizeof(double));
    outFile3.close();
    
    cout << "Arrays successfully saved to *.bin files (" << type << ")" << endl;
}

void clear_benchmark_data(double avrg_time[], double first_time[], double last_time[], int size) {
    fill_n(avrg_time, size, 0.0);
    fill_n(first_time, size, 0.0);
    fill_n(last_time, size, 0.0);
}
void print_search_info(int try_num, int x, int index, chrono::duration<double> elapsed) {
    cout << "  Try " << try_num << ": searching for " << x << endl;
    cout << "    Found at index: " << index << endl;
    cout << "    Time: " << fixed << setprecision(8) << elapsed.count() << " sec" << endl;
}
void bench_branchy_search(int n, int *index, int repeat = 1e6, int max_pow = 25) {
    chrono::duration<double> v, z, time_first, time_last;
    double avrg_time[25] = {0}; 
    double first_time[25] = {0};
    double last_time[25] = {0};
    
    cout << endl << "branchy_search" << endl;
    cout << "repeat each n = " << repeat << " times;" << endl;
    cout << "max power = 2^" << max_pow << endl << endl;
    
    srand(0);
    for (int i = 1; i <= max_pow; i++) {
        n = pow2(i);
        int* arr = new int[n];
        fillmas(arr, n);
        sorted_array<int, int> sa(arr, n);
        
        // cout << "Array size: 2^" << i << " (" << n << ")" << endl;
        // sa.print_array();
        
        z = std::chrono::seconds(0);
        for (int try1 = 0; try1 < repeat; try1++) {
            int x = rand() % n;
            auto start = chrono::steady_clock::now();
            *index = sa.branchy_search(x);
            auto end = chrono::steady_clock::now();
            chrono::duration<double> elapsed = end - start;
            
            if(try1 == 0) time_first = elapsed;
            if(try1 == repeat - 1) time_last = elapsed;
            z += elapsed;
            // print_search_info(try1 + 1, x, *index, elapsed);
        }
        
        auto result = z / repeat;
        cout << "n = 2^" << i << " (" << n << ")" << endl;
        cout << fixed << setprecision(8);
        cout << "  Total time: " << z.count() << " sec" << endl;
        cout << "  Average time: " << result.count() << " sec" << endl;
        cout << "  First try: " << time_first.count() << " sec" << endl;
        cout << "  Last try: " << time_last.count() << " sec" << endl << endl;
        
        avrg_time[i - 1] = result.count();
        first_time[i - 1] = time_first.count();
        last_time[i - 1] = time_last.count();
        
        delete[] arr;
    }
    create_bin_files(avrg_time, first_time, last_time, max_pow, "branchyBin");
    clear_benchmark_data(avrg_time, first_time, last_time, max_pow);
}

void bench_branchfree_search(int n, int *index, int repeat = 1e6, int max_pow = 25) {
    chrono::duration<double> v, z, time_first, time_last;
    double avrg_time[25] = {0}; 
    double first_time[25] = {0};
    double last_time[25] = {0};
    
    cout << endl << "branchfree_search" << endl;
    cout << "repeat each n = " << repeat << " times;" << endl;
    cout << "max power = 2^" << max_pow << endl << endl;
    
    srand(0);
    for (int i = 1; i <= max_pow; i++) {
        n = pow2(i);
        int* arr = new int[n];
        fillmas(arr, n);
        sorted_array<int, int> sa(arr, n);
        // cout << "Array size: 2^" << i << " (" << n << ")" << endl;
        // sa.print_array();
        z = std::chrono::seconds(0);
        for (int try1 = 0; try1 < repeat; try1++) {
            int x = rand() % n;
            auto start = chrono::steady_clock::now();
            *index = sa._branchfree_search(x);
            auto end = chrono::steady_clock::now();
            chrono::duration<double> elapsed = end - start;
            
            if(try1 == 0) time_first = elapsed;
            if(try1 == repeat - 1) time_last = elapsed;
            z += elapsed;
            // print_search_info(try1 + 1, x, *index, elapsed);
        }
        
        auto result = z / repeat;
        cout << "n = 2^" << i << " (" << n << ")" << endl;
        cout << fixed << setprecision(8);
        cout << "  Total time: " << z.count() << " sec" << endl;
        cout << "  Average time: " << result.count() << " sec" << endl;
        cout << "  First try: " << time_first.count() << " sec" << endl;
        cout << "  Last try: " << time_last.count() << " sec" << endl << endl;
        
        avrg_time[i - 1] = result.count();
        first_time[i - 1] = time_first.count();
        last_time[i - 1] = time_last.count();
        
        delete[] arr;
    }
    create_bin_files(avrg_time, first_time, last_time, max_pow, "branchfreeBin");
    clear_benchmark_data(avrg_time, first_time, last_time, max_pow);
}
void bench_branchfreePref2_search(int n, int *index, int repeat = 1e6, int max_pow = 25) {
    chrono::duration<double> v, z, time_first, time_last;
    double avrg_time[25] = {0}; 
    double first_time[25] = {0};
    double last_time[25] = {0};
    
    cout << endl << "branchfree_search" << endl;
    cout << "repeat each n = " << repeat << " times;" << endl;
    cout << "max power = 2^" << max_pow << endl << endl;
    
    srand(0);
    for (int i = 1; i <= max_pow; i++) {
        n = pow2(i);
        int* arr = new int[n];
        fillmas(arr, n);
        sorted_array<int, int> sa(arr, n);
        // cout << "Array size: 2^" << i << " (" << n << ")" << endl;
        // sa.print_array();
        z = std::chrono::seconds(0);
        for (int try1 = 0; try1 < repeat; try1++) {
            int x = rand() % n;
            auto start = chrono::steady_clock::now();
            *index = sa._branchfreepref_search(x);
            auto end = chrono::steady_clock::now();
            chrono::duration<double> elapsed = end - start;
            
            if(try1 == 0) time_first = elapsed;
            if(try1 == repeat - 1) time_last = elapsed;
            z += elapsed;
            // print_search_info(try1 + 1, x, *index, elapsed);
        }
        
        auto result = z / repeat;
        cout << "n = 2^" << i << " (" << n << ")" << endl;
        cout << fixed << setprecision(8);
        cout << "  Total time: " << z.count() << " sec" << endl;
        cout << "  Average time: " << result.count() << " sec" << endl;
        cout << "  First try: " << time_first.count() << " sec" << endl;
        cout << "  Last try: " << time_last.count() << " sec" << endl << endl;
        
        avrg_time[i - 1] = result.count();
        first_time[i - 1] = time_first.count();
        last_time[i - 1] = time_last.count();
        
        delete[] arr;
    }
    create_bin_files(avrg_time, first_time, last_time, max_pow, "branchfreePref2Bin");
    clear_benchmark_data(avrg_time, first_time, last_time, max_pow);
}

void bench_eytzinger_branchy_search(int n, int *index, int repeat = 1e6, int max_pow = 25) {
    chrono::duration<double> v, z, time_first, time_last;
    double avrg_time[25] = {0}; 
    double first_time[25] = {0};
    double last_time[25] = {0};
    
    cout << endl << "eytzinger_branchy_search" << endl;
    cout << "repeat each n = " << repeat << " times;" << endl;
    cout << "max power = 2^" << max_pow << endl << endl;
    
    srand(0);
    for (int i = 1; i <= max_pow; i++) {
        n = pow2(i);
        vector<int> arr;
        fillvec(arr, n);
        auto cbstRoot = BuildCompleteBST(arr, 0, arr.size() - 1);
        vector<int> levelOrder = LevelOrderToArray(cbstRoot);
        
        int* levelOrderArray = new int[levelOrder.size()];
        copy(levelOrder.begin(), levelOrder.end(), levelOrderArray);
        
        sorted_array<int, int> sa(levelOrderArray, levelOrder.size());
        // cout << "Array size: 2^" << i << " (" << n << ")" << endl;
        // sa.print_array();
        z = std::chrono::seconds(0);
        for (int try1 = 0; try1 < repeat; try1++) {
            int x = rand() % n;
            auto start = chrono::steady_clock::now();
            *index = sa.eytz_branchy_search(x);
            auto end = chrono::steady_clock::now();
            chrono::duration<double> elapsed = end - start;
            
            if(try1 == 0) time_first = elapsed;
            if(try1 == repeat - 1) time_last = elapsed;
            z += elapsed;
            // print_search_info(try1 + 1, x, *index, elapsed);
        }
        
        auto result = z / repeat;
        cout << "n = 2^" << i << " (" << n << ")" << endl;
        cout << fixed << setprecision(8);
        cout << "  Total time: " << z.count() << " sec" << endl;
        cout << "  Average time: " << result.count() << " sec" << endl;
        cout << "  First try: " << time_first.count() << " sec" << endl;
        cout << "  Last try: " << time_last.count() << " sec" << endl << endl;
        
        avrg_time[i - 1] = result.count();
        first_time[i - 1] = time_first.count();
        last_time[i - 1] = time_last.count();
        
        delete[] levelOrderArray;
        arr.clear();
        levelOrder.clear();
    }
    create_bin_files(avrg_time, first_time, last_time, max_pow, "eytzingerBin");
    clear_benchmark_data(avrg_time, first_time, last_time, max_pow);
}

void bench_eytzinger_branchfree_search(int n, int *index, int repeat = 1e6, int max_pow = 25) {
    chrono::duration<double> v, z, time_first, time_last;
    double avrg_time[25] = {0}; 
    double first_time[25] = {0};
    double last_time[25] = {0};
    
    cout << endl << "eytzinger_branchfree_search" << endl;
    cout << "repeat each n = " << repeat << " times;" << endl;
    cout << "max power = 2^" << max_pow << endl << endl;
    
    srand(0);
    for (int i = 1; i <= max_pow; i++) {
        n = pow2(i);
        vector<int> arr;
        fillvec(arr, n);
        auto cbstRoot = BuildCompleteBST(arr, 0, arr.size() - 1);
        vector<int> levelOrder = LevelOrderToArray(cbstRoot);
        
        int* levelOrderArray = new int[levelOrder.size()];
        copy(levelOrder.begin(), levelOrder.end(), levelOrderArray);
        
        sorted_array<int, int> sa(levelOrderArray, levelOrder.size());
        // cout << "Array size: 2^" << i << " (" << n << ")" << endl;
        // sa.print_array();
        z = std::chrono::seconds(0);
        for (int try1 = 0; try1 < repeat; try1++) {
            int x = rand() % n;
            auto start = chrono::steady_clock::now();
            *index = sa.eytz_branchfree_search(x);
            auto end = chrono::steady_clock::now();
            chrono::duration<double> elapsed = end - start;
            
            if(try1 == 0) time_first = elapsed;
            if(try1 == repeat - 1) time_last = elapsed;
            z += elapsed;
            // print_search_info(try1 + 1, x, *index, elapsed);
        }
        
        auto result = z / repeat;
        cout << "n = 2^" << i << " (" << n << ")" << endl;
        cout << fixed << setprecision(8);
        cout << "  Total time: " << z.count() << " sec" << endl;
        cout << "  Average time: " << result.count() << " sec" << endl;
        cout << "  First try: " << time_first.count() << " sec" << endl;
        cout << "  Last try: " << time_last.count() << " sec" << endl << endl;
        
        avrg_time[i - 1] = result.count();
        first_time[i - 1] = time_first.count();
        last_time[i - 1] = time_last.count();
        
        delete[] levelOrderArray;
        arr.clear();
        levelOrder.clear();
    }
    create_bin_files(avrg_time, first_time, last_time, max_pow, "eytzingerFreeBin");
    clear_benchmark_data(avrg_time, first_time, last_time, max_pow);
}


int main() {
    int index;
    int n;
    
     bench_branchy_search(n, &index);
     bench_branchfree_search(n, &index);
     bench_branchfreePref2_search(n, &index);
     bench_eytzinger_branchy_search(n, &index);
     bench_eytzinger_branchfree_search(n, &index);
    
    return 0;
}
