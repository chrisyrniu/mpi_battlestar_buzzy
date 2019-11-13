/*
Author: <Yaru Niu>
Class: ECE6122
Last Date Modified: <2019-11-13>
Description: Solution to Homework 5
*/

#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <math.h>

#include "mpi.h"

using namespace std;

// Main
int main(int argc, char** argv){

    //***** Get information from the input file ******//
    // Define variables to read the file
    fstream input;
    double number;
    int flag1 = 0;
    int flag2 = 0;
    int timeLimit;
    double maxThrust;
    int rows = 8;
    int columns = 7;
    double matrix[8][7];

    // Read the file, and store the maximum time, maximum thrust
    // and initial positions and velocities of Buzzy and Yellow Jackets
    const char* filepath = "in.dat";
    input.open(filepath, ios::in);
    if (!input){
        cout << "Unable to open file";
    }
    while(input >> number){
        if (flag1 == 0){
            timeLimit = number;
            flag1 = flag1 + 1;
        }
        else if (flag1 == 1){
            maxThrust = number;
            flag1 = flag1 + 1;
        }
        else{
            matrix[flag2/columns][flag2%columns] = number;
            flag2 = flag2 + 1;
        }
    }
    input.close();

    //***************** MPI ********************//

    int  numtasks, rank, rc;

    // Initialize the MPI execution environment
    rc = MPI_Init(&argc, &argv);

    // Check if initialize successfully, otherwise abort
    if (rc != MPI_SUCCESS)
    {
        printf("Error starting MPI program. Terminating.\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }

    // Return the total number of MPI processes in the specified communicator (MPI_COMM_WORLD)
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    // Return the rank of the calling MPI process within the specified communicator (MPI_COMM_WORLD)
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Get the initial velocity and direction
    double vel = matrix[rank][3];
    double decX = matrix[rank][4];
    double decY = matrix[rank][5];
    double decZ = matrix[rank][6];
    // Get the initial position
    double xPos = matrix[rank][0];
    double yPos = matrix[rank][1];
    double zPos = matrix[rank][2]; 

    // Initialize msg (1D array) which stores information of Buzzy or Yellow Jackets
    double msg[8];
    msg[0] = rank;
    msg[1] = 1;
    msg[2] = xPos;
    msg[3] = yPos;
    msg[4] = zPos;
    msg[5] = 0;
    msg[6] = 0;
    msg[7] = 0;
    // Define msgRec (1D array) which is used for storing received information
    double msgRec[8];
    // Define msgCo (2D array, composite msg) which stores all information to output
    double msgCo[8][8];
    // Define new MPI datatype rowType (1D array)
    MPI_Datatype rowType;
    MPI_Type_contiguous(8, MPI_DOUBLE, &rowType);
    MPI_Type_commit(&rowType);
    // Define new MPI datatype matrixType (2D array)
    MPI_Datatype matrixType;
    MPI_Type_contiguous(64, MPI_DOUBLE, &matrixType);
    MPI_Type_commit(&matrixType);

    // Loops within the time limit
    for (int round = 0; round < timeLimit; ++round)
    {
        // Update position using the velocity and direction
        msg[2] = msg[2] + vel*decX;
        msg[3] = msg[3] + vel*decY;
        msg[4] = msg[4] + vel*decZ;

        // Main process
        if (rank == 0)
        {
            // Update information of Buzzy
            for (int m = 0; m < 8; ++m)
            {
                msgCo[0][m] = msg[m];
            }

            // Required variable for receive routines
            MPI_Status status;
            // Receive messages (1D array) from Yellow Jackets, get composite msg
            for (int i = 1; i < 8; ++i)
            {
                rc = MPI_Recv(msgRec, 8, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
                for (int j = 0; j < 8; ++j)
                {
                    msgCo[i][j] = msgRec[j];
                }
            }
            // Send the composite msg to all the Yellow Jackets
            for (int k = 1; k < 8; ++k)
            {
                rc = MPI_Send(msgCo, 1, matrixType, k, 0, MPI_COMM_WORLD);
            }
            // Output the information (rankID, status, location and thrust) of each Yellow Jackets
            for (int a = 1; a < 8; ++a)
            {
                for (int b = 0; b < 8; ++b)
                {
                    if (b < 2)
                    {
                        cout << int(msgCo[a][b]) << ", "; 
                    }
                    else if (b < 7)
                    {
                        cout << scientific << setprecision(6) << msgCo[a][b] << ", "; 
                    }
                    else
                        cout << scientific << setprecision(6) << msgCo[a][b] << endl;
                }
            }
        }
        // Yellow Jacket process
        else
        {
            // Send the message to Buzzy
            rc = MPI_Send(&msg[0], 1, rowType, 0, 0, MPI_COMM_WORLD);         
            // Required variable for receive routines
            MPI_Status status;
            // Receive composite msg from Buzzy
            rc = MPI_Recv(msgCo, 64, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        }        
    }

    // Free datatype when done using them
    MPI_Type_free(&rowType);
    MPI_Type_free(&matrixType);
    // Finalize
    MPI_Finalize();
}
