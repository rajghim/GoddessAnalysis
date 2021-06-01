#define Analysis_cxx
#include "Analysis.h"
#include "EnergyLoss.h"
#include <TRandom3.h>
#include <string>
#include "EffectiveThickness.h"
#include "json/json.h"


TChain* MakeChain();

int main() {
    TChain *chain = MakeChain();
    Analysis t(chain);
    t.Loop();
}

TChain* MakeChain() {
    auto *chain = new TChain("data");

    //Read and Parse config.json
    Json::Value config;
    std::ifstream config_stream("config.json");
    ASSERT_WITH_MESSAGE(config_stream.is_open(),
                        "Could not find 'config.json'\n");
    config_stream >> config;
    config_stream.close();

    TString InputPath = config["InputPath"].asString();
    TString InputFilePrefix = config["InputFilePrefix"].asString();

    for(int i = 0; i < config["Runs"].size(); i++) {
        TString RunNumber = config["Runs"][i].asString();
        //Add all runs to the chain
        chain->Add(InputPath + InputFilePrefix + RunNumber + "_combined.root"); 
    } 
    return chain;
}
