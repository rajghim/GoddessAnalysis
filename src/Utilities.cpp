#include "Utilities.h"
#include <iostream>
#include <fstream>


/*const std::string red("\033[1;31m");
const std::string green("\033[1;32m]");
const std::string yellow("\033[1;33m]");
const std::string blue("\033[1;34m]");
const std::string magenta("\033[1;35m]");
const std::string cyan("\033[1;36m]");
const std::string reset("\033[0m");*/
//This function outputs effective thickness
void DownstreamWelcomeMessage(){ 

std::cout << green << "               ___         _    _               _             _         _" << std::endl;                                                           
std::cout << "              / __|___  __| |__| |___ ______   /_\\  _ _  __ _| |_  _ __(_)___" << std::endl;                                                      
std::cout << "             | (_ / _ \\/ _` / _` / -_|_-<_-<  / _ \\| ' \\/ _` | | || (_-< (_-<" << std::endl;                                                      
std::cout << "              \\___\\___/\\__,_\\__,_\\___/__/__/ /_/ \\_\\_||_\\__,_|_|\\_, /__/_/__/" << std::endl;                                                       
std::cout << "              ___                     _                        _|__/                 _    _______  ______    _     ___ ___ _  ____ " << std::endl; 
std::cout << "             |   \\ _____ __ ___ _  __| |_ _ _ ___ __ _ _ __   | _ ) __ _ _ _ _ _ ___| |  / / __\\ \\/ /__ /  _| |_  | _ ) _ ) |/  \\ \\ " << std::endl;
std::cout << "             | |) / _ \\ V  V / ' \\(_-<  _| '_/ -_) _` | '  \\  | _ \\/ _` | '_| '_/ -_) | | |\\__ \\>  < |_ \\ |_   _| | _ \\ _ \\ | () | |" << std::endl;
std::cout << "             |___/\\___/\\_/\\_/|_||_/__/\\__|_| \\___\\__,_|_|_|_| |___/\\__,_|_| |_| \\___|_| | ||___/_/\\_\\___/   |_|   |___/___/_|\\__/| |" << std::endl;
std::cout << "                                                                                         \\_\\                                    /_/ " << reset << std::endl;   

}    


void UpstreamQQQ5WelcomeMessage(){

std::cout << green << "                       ___         _    _               _             _         _" << std::endl;     
std::cout << "                      / __|___  __| |__| |___ ______   /_\\  _ _  __ _| |_  _ __(_)___" << std::endl; 
std::cout << "                     | (_ / _ \\/ _` / _` / -_|_-<_-<  / _ \\| ' \\/ _` | | || (_-< (_-<" << std::endl; 
std::cout << "                      \\___\\___/\\__,_\\__,_\\___/__/__/ /_/ \\_\\_||_\\__,_|_|\\_, /__/_/__/" << std::endl; 
std::cout << "                     | | | |_ __ __| |_ _ _ ___ __ _ _ __    __ _ __ _ _|__/ __|" << std::endl;      
std::cout << "                     | |_| | '_ (_-<  _| '_/ -_) _` | '  \\  / _` / _` / _` |__ \\" << std::endl;      
std::cout << "                      \\___/| .__/__/\\__|_| \\___\\__,_|_|_|_| \\__, \\__, \\__, |___/ " << std::endl;     
std::cout << "                           |_|                                 |_|  |_|  |_|  " << reset << std::endl; 



}


void UpstreamSX3WelcomeMessage(){

std::cout << green << "                       ___         _    _               _             _         _ " << std::endl;   
std::cout << "                      / __|___  __| |__| |___ ______   /_\\  _ _  __ _| |_  _ __(_)___" << std::endl;
std::cout << "                     | (_ / _ \\/ _` / _` / -_|_-<_-<  / _ \\| ' \\/ _` | | || (_-< (_-<" << std::endl;
std::cout << "                      \\___\\___/\\__,_\\__,_\\___/__/__/ /_/ \\_\\_||_\\__,_|_|\\_, /__/_/__/" << std::endl;
std::cout << "                      _   _         _                        _____  ____|__/" << std::endl;         
std::cout << "                     | | | |_ __ __| |_ _ _ ___ __ _ _ __   / __\\ \\/ /__ / " << std::endl;          
std::cout << "                     | |_| | '_ (_-<  _| '_/ -_) _` | '  \\  \\__ \\>  < |_ \\ " << std::endl;          
std::cout << "                      \\___/| .__/__/\\__|_| \\___\\__,_|_|_|_| |___/_/\\_\\___/ " << std::endl;          
std::cout << "                           |_|                                             " << reset << std::endl;          

}       