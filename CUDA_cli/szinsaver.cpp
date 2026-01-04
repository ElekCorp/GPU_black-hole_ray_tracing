#include "szinsaver.h"
#include <iostream>
#include <fstream>

int datasaver(int8_t* szin,int SZELES,int MAGAS, char* filename)
{
    std::ofstream wf(filename, std::ios::out | std::ios::binary);
    if (!wf) {
        std::cout << "Cannot open file!" << std::endl;
        return 1;
    }

    int m = MAGAS;
    int sz = SZELES;

    wf.write((char*)&m, sizeof(int));
    wf.write((char*)&sz, sizeof(int));

    for (int i = 0; i < MAGAS * SZELES; i++)
    {
        wf.write((char*)&szin[i], sizeof(int8_t));
    }
    wf.close();
    if (!wf.good())
    {
        std::cout << "Error occurred at writing time!" << std::endl;
        return 1;
    }
    return 0;

}


int datasaver(int8_t* szin, int SZELES, int MAGAS, std::string& filename)
{
    std::ofstream wf(filename, std::ios::out | std::ios::binary);
    if (!wf) {
        std::cout << "Cannot open file!" << std::endl;
        return 1;
    }

    int m = MAGAS;
    int sz = SZELES;

    wf.write((char*)&m, sizeof(int));
    wf.write((char*)&sz, sizeof(int));

    for (int i = 0; i < MAGAS * SZELES; i++)
    {
        wf.write((char*)&szin[i], sizeof(int8_t));
    }
    wf.close();
    if (!wf.good())
    {
        std::cout << "Error occurred at writing time!" << std::endl;
        return 1;
    }
    return 0;

}

