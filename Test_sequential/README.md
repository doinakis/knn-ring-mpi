# knnring
Implementation of knn (V0. sequential)  
## Run the code with its own main function  
**In order to run the code with the provided main function you need to go to the knnring directory through the bash and do the following**:
1.  change directory and go to the **cd path/knn-ring-mpi/Test_sequential/knnring** where path is where you saved the repository in your own computer 
2.  type **make** to make the .a files and the executable main function i included
3.  cd to the src folder wich is located in the knnring folder and type **cd path/knn-ring-mpi/Test_sequential/knnring/src**
4. type ./main_sequential to run the main function and it will print how much time it took to run in your computer.
## Run the code with the tester 
To run the code with the tester :
1. repeat the first two steps of the previous instructions 
2. cd to the Test_sequential folder again and make a tar.gz file of knnring using **tar -czvf code.tar.gz knnring/**
3. type **make** in the same folder and the tester will run and print the validation results. 
### Important 
I already provide you the tester that is uploaded at e-learning but if you want you can include yours (if you dont trust me :P) just remember to **delete the test_mpi part to test only the sequential part.** 
