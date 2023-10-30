#ifndef __FILE_CTRL_H__
#define __FILE_CTRL_H__

#include <iostream>
#include <fstream>
#include <cstdlib>
#include "serialize.h"

using namespace std;

void split_filename(const string& filepath, string &dir, string &filename)
{   
    size_t found = filepath.find_last_of("/");
    if (found ==  string::npos)
    {
        dir = "";
        filename = filepath;
    }
    else
    {
        dir = filepath.substr(0, found);
        filename = filepath.substr(found + 1);
    }

    cout << " path: " << dir << '\n';
    cout << " file: " << filename << '\n';
}

int make_dir(const string& filepath)
{
    string dir;
    string filename;
    split_filename(filepath, dir, filename);
    if (dir == "")
    {
        return 0;
    }

    string cmd = "mkdir -p " + dir;
    return system(cmd.c_str());
}

/// Save a serialized file
template <class T>
static void save_file(const string filename, const T& output)
{
    ofstream outfile(filename, ios::binary);

    if (!outfile.eof() && !outfile.fail())
    {
        StreamType res;
        serialize(output, res);
        outfile.write(reinterpret_cast<char*>(&res[0]), res.size());
        outfile.close();
        res.clear();
        cout << "Save file successfully: " << filename << '\n';
    }
    else
    {
        cout << "Save file failed: " + filename << '\n';
        exit(1);
    }
}

/// Load a serialized file
template <class T>
static void load_file(const string filename, T& input)
{
    ifstream infile(filename, ios::binary);

    if (!infile.eof() && !infile.fail())
    {
        infile.seekg(0, ios_base::end);
        const streampos fileSize = infile.tellg();
        infile.seekg(0, ios_base::beg);
        vector<uint8_t> res(fileSize);
        infile.read(reinterpret_cast<char*>(&res[0]), fileSize);
        infile.close();
        input.clear();
        auto it = res.cbegin();
        input = deserialize<T>(it, res.cend());
        res.clear();
    }
    else
    {
        cout << "Cannot open file: " + filename << '\n';
        exit(1);
    }
}

/// Save graph structure to a file
template <class T>
void save_serialized_graph(const string file_name, const T& graph)
{
    make_dir(file_name);
    save_file(file_name, graph);
}

/// Load graph structure from a file
template <class T>
void load_serialized_graph(const string file_name, T& graph)
{
    load_file(file_name, graph);
}

bool seraizlied_graph_exist(const string file_name)
{
    ifstream infile(file_name, ios::binary);

    if (!infile.is_open())
    {
        return false;
    }

    infile.close();
    return true;
}
#endif