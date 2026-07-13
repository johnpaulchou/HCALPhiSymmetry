#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/DetId/interface/DetId.h"

#include <iomanip>
#include <iostream>
#include <vector>



class PrintHFGeometry : public edm::one::EDAnalyzer<> {
public:
  explicit PrintHFGeometry(edm::ParameterSet const&)
      : caloGeomToken_(esConsumes<CaloGeometry, CaloGeometryRecord>()) {}

  void analyze(edm::Event const&, edm::EventSetup const& iSetup) override {
    auto const& caloGeometry = iSetup.getData(caloGeomToken_);

    auto const& detIds =
        caloGeometry.getValidDetIds(DetId::Hcal, HcalForward);

    std::cout << "{\n";
    std::cout << "  \"hf_cells\": [\n";

    bool first = true;
    
    for (auto const& id : detIds) {
      HcalDetId hid(id);

      auto const cell = caloGeometry.getGeometry(id);
      if (!cell.get()) {
        std::cout << "No CaloCellGeometry for rawId=" << id.rawId() << "\n";
        continue;
      }

      GlobalPoint const pos = cell->getPosition();

      double const frontArea = quadArea(cell->getCorners(), 0, 1, 2, 3);
      double const backArea  = quadArea(cell->getCorners(), 4, 5, 6, 7);
      double const avgArea   = 0.5 * (frontArea + backArea);

      if (!first) {
        std::cout << ",\n";
      }
      first = false;

      std::cout << "    {\n";
      std::cout << "      \"rawId\": " << id.rawId() << ",\n";
      std::cout << "      \"ieta\": " << hid.ieta() << ",\n";
      std::cout << "      \"iphi\": " << hid.iphi() << ",\n";
      std::cout << "      \"depth\": " << hid.depth() << ",\n";
      std::cout << "      \"eta\": " << pos.eta() << ",\n";
      std::cout << "      \"phi\": " << pos.phi() << ",\n";
      std::cout << "      \"x_cm\": " << pos.x() << ",\n";
      std::cout << "      \"y_cm\": " << pos.y() << ",\n";
      std::cout << "      \"z_cm\": " << pos.z() << ",\n";
      std::cout << "      \"frontArea_cm2\": " << frontArea << ",\n";
      std::cout << "      \"backArea_cm2\": " << backArea << ",\n";
      std::cout << "      \"avgArea_cm2\": " << avgArea << ",\n";
      std::cout << "      \"etaSpan\": " << cell->etaSpan() << ",\n";
      std::cout << "      \"phiSpan\": " << cell->phiSpan() << "\n";
      std::cout << "    }";
    }

    std::cout << "\n";
    std::cout << "  ]\n";
    std::cout << "}\n";
  }

private:
  using Corners = CaloCellGeometry::CornersVec;

  static GlobalVector diff(GlobalPoint const& a, GlobalPoint const& b) {
    return GlobalVector(a.x() - b.x(), a.y() - b.y(), a.z() - b.z());
  }

  static double triangleArea(GlobalPoint const& a,
                             GlobalPoint const& b,
                             GlobalPoint const& c) {
    GlobalVector const ab = diff(b, a);
    GlobalVector const ac = diff(c, a);
    return 0.5 * ab.cross(ac).mag();
  }

  static double quadArea(Corners const& c,
                         unsigned int i0,
                         unsigned int i1,
                         unsigned int i2,
                         unsigned int i3) {
    // Split quadrilateral into two triangles.
    return triangleArea(c[i0], c[i1], c[i2]) +
           triangleArea(c[i0], c[i2], c[i3]);
  }

  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;
};

DEFINE_FWK_MODULE(PrintHFGeometry);
