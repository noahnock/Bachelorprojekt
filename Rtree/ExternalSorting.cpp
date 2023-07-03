// C++ program to implement
// external sorting using
// merge sort
// by https://www.geeksforgeeks.org/external-sorting/ 02.07.2023
#include <fstream>
#include <utility>
#include <vector>
#include <iostream>
#include <cmath>
#include "./Rtree.h"
#include <float.h>
using namespace std;

auto sortRuleLambdaX = [] (rTreeValue b1, rTreeValue b2) -> bool
{
    double center1 = (b1.first.min_corner().get<0>() + b1.first.max_corner().get<0>()) / 2;
    double center2 = (b2.first.min_corner().get<0>() + b2.first.max_corner().get<0>()) / 2;
    return center1 < center2;
};

auto sortRuleLambdaY = [](rTreeValue b1, rTreeValue b2) -> bool {
    double center1 = (b1.first.min_corner().get<1>() + b1.first.max_corner().get<1>()) / 2;
    double center2 = (b2.first.min_corner().get<1>() + b2.first.max_corner().get<1>()) / 2;
    return center1 < center2;
};

struct MinHeapNode {

    // The element to be stored
    rTreeValue element;

    // index of the array from which
    // the element is taken
    int i;
};

// Prototype of a utility function
// to swap two min heap nodes
void swap(MinHeapNode* x, MinHeapNode* y);

// A class for Min Heap
class MinHeap {

    // pointer to array of elements in heap
    MinHeapNode* harr;

    // size of min heap
    int heap_size;

    // dimension of the sorting (0 is x, 1 is y)
    size_t dim;

public:

    // Constructor: creates a min
    // heap of given size
    MinHeap(MinHeapNode a[], int size, size_t dim);

    // to heapify a subtree with
    // root at given index
    void MinHeapify(int);

    // to get index of left child
    // of node at index i
    static int left(int i) { return (2 * i + 1); }

    // to get index of right child
    // of node at index i
    static int right(int i) { return (2 * i + 2); }

    // to get the root
    MinHeapNode getMin() { return harr[0]; }

    // to replace root with new node
    // x and heapify() new root
    void replaceMin(MinHeapNode x)
    {
        harr[0] = x;
        MinHeapify(0);
    }
};

// Constructor: Builds a heap from
// a given array a[] of given size
MinHeap::MinHeap(MinHeapNode a[], int size, size_t dim)
{
    heap_size = size;
    harr = a; // store address of array
    this->dim = dim;
    int i = (heap_size - 1) / 2;
    while (i >= 0) {
        MinHeapify(i);
        i--;
    }
}

// A recursive method to heapify
// a subtree with root
// at given index. This method
// assumes that the
// subtrees are already heapified
void MinHeap::MinHeapify(int i)
{
    int l = left(i);
    int r = right(i);
    int smallest = i;

    if (dim == 0) {
        if (l < heap_size && sortRuleLambdaX(harr[l].element, harr[i].element))
            smallest = l;

        if (r < heap_size
            && sortRuleLambdaX(harr[r].element, harr[smallest].element))
            smallest = r;
    } else {
        if (l < heap_size && sortRuleLambdaY(harr[l].element, harr[i].element))
            smallest = l;

        if (r < heap_size
            && sortRuleLambdaY(harr[r].element, harr[smallest].element))
            smallest = r;
    }

    if (smallest != i) {
        swap(&harr[i], &harr[smallest]);
        MinHeapify(smallest);
    }
}

// A utility function to swap two elements
void swap(MinHeapNode* x, MinHeapNode* y)
{
    MinHeapNode temp = *x;
    *x = *y;
    *y = temp;
}

// Merges two subarrays of arr[].
// First subarray is arr[l..m]
// Second subarray is arr[m+1..r]
void merge(rTreeValue arr[], long long l, long long m, long long r, size_t dim)
{
    long long i, j, k;
    long long n1 = m - l + 1;
    long long n2 = r - m;

    /* create temp arrays */
    rTreeValue L[n1], R[n2];

    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++)
        L[i] = arr[l + i];
    for (j = 0; j < n2; j++)
        R[j] = arr[m + 1 + j];

    /* Merge the temp arrays back into arr[l..r]*/
    // Initial index of first subarray
    i = 0;

    // Initial index of second subarray
    j = 0;

    // Initial index of merged subarray
    k = l;
    while (i < n1 && j < n2) {
        if (dim == 0) {
            if (sortRuleLambdaX(L[i], R[j]))
                arr[k++] = L[i++];
            else
                arr[k++] = R[j++];
        } else {
            if (sortRuleLambdaY(L[i], R[j]))
                arr[k++] = L[i++];
            else
                arr[k++] = R[j++];
        }
    }

    /* Copy the remaining elements of L[],
        if there are any */
    while (i < n1)
        arr[k++] = L[i++];

    /* Copy the remaining elements of R[],
        if there are any */
    while (j < n2)
        arr[k++] = R[j++];
}

/* l is for left index and r is right index of the
sub-array of arr to be sorted */
void mergeSort(rTreeValue arr[], long long l, long long r, size_t dim)
{
    if (l < r) {
        // Same as (l+r)/2, but avoids overflow for
        // large l and h
        long long m = l + (r - l) / 2;

        // Sort first and second halves
        mergeSort(arr, l, m, dim);
        mergeSort(arr, m + 1, r, dim);

        merge(arr, l, m, r, dim);
    }
}

// Merges k sorted files. Names of files are assumed
// to be 1, 2, 3, ... k
void mergeFiles(const std::string& output_file, int k, size_t dim)
{
    std::vector<ifstream> in(k);
    for (int i = 0; i < k; i++) {

        // convert i to string
        std::string fileName = "external-sorting-tmp-" + to_string(i);

        // Open output files in read mode.
        in[i] = ifstream(fileName, std::ios::binary);
    }

    // FINAL OUTPUT FILE
    ofstream out = ofstream(output_file, std::ios::binary);

    // Create a min heap with k heap
    // nodes. Every heap node
    // has first element of scratch
    // output file
    auto* harr = new MinHeapNode[k];
    int i;
    for (i = 0; i < k; i++) {
        // break if no output file is empty and
        // index i will be no. of input files
        long long pos = in[i].tellg();
        in[i].seekg (0, std::ifstream::end);
        long long fileLength = in[i].tellg();
        in[i].seekg (pos, std::ifstream::beg);
        if (pos < fileLength) {
            double minX;
            double minY;
            double maxX;
            double maxY;
            uint64_t id;
            in[i].read(reinterpret_cast<char*>(&minX), sizeof(double));
            in[i].read(reinterpret_cast<char*>(&minY), sizeof(double));
            in[i].read(reinterpret_cast<char*>(&maxX), sizeof(double));
            in[i].read(reinterpret_cast<char*>(&maxY), sizeof(double));
            in[i].read(reinterpret_cast<char*>(&id), sizeof(uint64_t));

            auto box = make<boxGeo>(make<pointGeo>(minX, minY), make<pointGeo>(maxX, maxY));;
            rTreeValue boxWithId = std::make_pair(box, id);
            harr[i].element = boxWithId;
        } else {
            break;
        }

        // Index of scratch output file
        harr[i].i = i;
    }
    // Create the heap
    MinHeap hp(harr, i, dim);

    int count = 0;

    // Now one by one get the
    // minimum element from min
    // heap and replace it with
    // next element.
    // run till all filled input
    // files reach EOF
    while (count != i) {
        // Get the minimum element
        // and store it in output file
        MinHeapNode root = hp.getMin();

        double minXWrite = root.element.first.min_corner().get<0>();
        double minYWrite = root.element.first.min_corner().get<1>();
        double maxXWrite = root.element.first.max_corner().get<0>();
        double maxYWrite = root.element.first.max_corner().get<1>();

        out.write(reinterpret_cast<const char *>(&minXWrite), sizeof(double));
        out.write(reinterpret_cast<const char *>(&minYWrite), sizeof(double));
        out.write(reinterpret_cast<const char *>(&maxXWrite), sizeof(double));
        out.write(reinterpret_cast<const char *>(&maxYWrite), sizeof(double));
        out.write(reinterpret_cast<const char *>(&root.element.second), sizeof(long long));

        // Find the next element that
        // will replace current
        // root of heap. The next element
        // belongs to same
        // input file as the current min element.
        long long pos = in[root.i].tellg();
        in[root.i].seekg (0, std::ifstream::end);
        long long fileLength = in[root.i].tellg();
        in[root.i].seekg (pos, std::ifstream::beg);
        if (pos < fileLength) {
            double minX;
            double minY;
            double maxX;
            double maxY;
            uint64_t id;
            in[root.i].read(reinterpret_cast<char*>(&minX), sizeof(double));
            in[root.i].read(reinterpret_cast<char*>(&minY), sizeof(double));
            in[root.i].read(reinterpret_cast<char*>(&maxX), sizeof(double));
            in[root.i].read(reinterpret_cast<char*>(&maxY), sizeof(double));
            in[root.i].read(reinterpret_cast<char*>(&id), sizeof(uint64_t));

            auto box = make<boxGeo>(make<pointGeo>(minX, minY), make<pointGeo>(maxX, maxY));
            rTreeValue boxWithId = std::make_pair(box, id);
            root.element = boxWithId;
        } else {
            auto box = make<boxGeo>(make<pointGeo>(DBL_MAX, DBL_MAX), make<pointGeo>(DBL_MAX, DBL_MAX));
            rTreeValue boxWithId = std::make_pair(box, INT_MAX);
            root.element = boxWithId;
            count++;
        }

        // Replace root with next
        // element of input file
        hp.replaceMin(root);
    }

    // close input and output files
    for (int j = 0; j < k; j++)
        in[j].close();

    out.close();
}

// Using a merge-sort algorithm,
// create the initial runs
// and divide them evenly among
// the output files
void createInitialRuns(const std::string& input_file, uintmax_t run_size,
                       int num_ways, size_t dim)
{
    // For big input file
    ifstream in = ifstream(input_file, std::ios::binary);

    // output scratch files
    std::vector<ofstream> out(num_ways);
    for (int i = 0; i < num_ways; i++) {
        // convert i to string
        std::string fileName = "external-sorting-tmp-" + to_string(i);

        // Open output files in write mode.
        out[i] = ofstream(fileName, std::ios::binary);
    }

    // allocate a dynamic array large enough
    // to accommodate runs of size run_size
    auto* arr = (rTreeValue*)malloc(run_size * sizeof(rTreeValue));

    bool more_input = true;
    int next_output_file = 0;

    uintmax_t i;

    in.seekg (0, std::ifstream::end);
    long long fileLength = in.tellg();
    in.seekg (0, std::ifstream::beg);

    while (more_input) {
        // write run_size elements
        // into arr from input file
        for (i = 0; i < run_size; i++) {
            if (in.tellg() < fileLength) {
                /*double minX;
                double minY;
                double maxX;
                double maxY;
                uint64_t id;
                in.read(reinterpret_cast<char*>(&minX), sizeof(double));
                in.read(reinterpret_cast<char*>(&minY), sizeof(double));
                in.read(reinterpret_cast<char*>(&maxX), sizeof(double));
                in.read(reinterpret_cast<char*>(&maxY), sizeof(double));
                in.read(reinterpret_cast<char*>(&id), sizeof(uint64_t));

                auto box = make<boxGeo>(make<pointGeo>(minX, minY), make<pointGeo>(maxX, maxY));;
                rTreeValue boxWithId = std::make_pair(box, id);
                arr[i] = boxWithId;*/
            } else {
                more_input = false;
                break;
            }
        }

        // sort array using merge sort
        //mergeSort(arr, 0, i - 1, dim);

        // write the records to the
        // appropriate scratch output file
        // can't assume that the loop
        // runs to run_size
        // since the last run's length
        // may be less than run_size
        /*for (int j = 0; j < i; j++) {
            double minX = arr[j].first.min_corner().get<0>();
            double minY = arr[j].first.min_corner().get<1>();
            double maxX = arr[j].first.max_corner().get<0>();
            double maxY = arr[j].first.max_corner().get<1>();

            out[next_output_file].write(reinterpret_cast<const char *>(&minX), sizeof(double));
            out[next_output_file].write(reinterpret_cast<const char *>(&minY), sizeof(double));
            out[next_output_file].write(reinterpret_cast<const char *>(&maxX), sizeof(double));
            out[next_output_file].write(reinterpret_cast<const char *>(&maxY), sizeof(double));
            out[next_output_file].write(reinterpret_cast<const char *>(&arr[j].second), sizeof(long long));
        }*/

        next_output_file++;
    }

    // close input and output files
    for (int j = 0; j < num_ways; j++)
        out[j].close();

    in.close();
}

// For sorting data stored on disk
void externalSort(const std::string& input_file, const std::string& output_file, uintmax_t maxBuildingRamUsage, size_t dim)
{
    int num_ways = std::ceil(((double) std::filesystem::file_size(input_file)) / ((double) maxBuildingRamUsage));
    // read the input file,
    // create the initial runs,
    // and assign the runs to
    // the scratch output files
    createInitialRuns(input_file, maxBuildingRamUsage, num_ways, dim);

    // Merge the runs using
    // the K-way merging
    //mergeFiles(output_file, num_ways, dim);

    for (int j = 0; j < num_ways; j++) {
        std::remove(("external-sorting-tmp-" + to_string(j)).c_str());
    }
}

// Driver code
void run()
{

    // The size of each partition
    int run_size = 1000;

    int total_size = 10000;

    char input_file[] = "../input.txt";
    char output_file[] = "../output.txt";

    std::ofstream in = std::ofstream(input_file, std::ios::binary);

    srand(time(nullptr));

    // generate input
    std::vector<int> input;
    for (int i = 0; i < total_size; i++) {
        int r = rand();
        in.write(reinterpret_cast<const char *>(&r), sizeof(int));
        input.push_back(r);
    }

    // No. of Partitions of input file.
    int num_ways = std::ceil(((double) total_size) / run_size);

    in.close();

    externalSort(input_file, output_file,
                 run_size, 0);

    std::sort(input.begin(), input.end());
    ifstream result = std::ifstream( "../output.txt", std::ios::binary);
    std::string line;

    int i = 0;
    bool correct = true;

    result.seekg (0, std::ifstream::end);
    long long fileLength = result.tellg();
    result.seekg (0, std::ifstream::beg);
    while(result.tellg() < fileLength)
    {
        int currentFromFile;
        result.read(reinterpret_cast<char*>(&currentFromFile), sizeof(int));
        int currentFromVector = input[i];
        if (currentFromFile != currentFromVector) {
            correct = false;
        }
        i++;
    }

    std::cout << (correct ? "True" : "False") << std::endl;

    std::remove("../input.txt");
    std::remove("../output.txt");
    for (int j = 0; j < 40; j++) {
        std::remove(to_string(j).c_str());
    }

    result.close();
}

