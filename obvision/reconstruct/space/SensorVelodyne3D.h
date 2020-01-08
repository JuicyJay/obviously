#ifndef OBVISION_RECONSTRUCT_SPACE_SENSORVELODYNE3D_H_
#define OBVISION_RECONSTRUCT_SPACE_SENSORVELODYNE3D_H_

#include "obvision/reconstruct/Sensor.h"
#include "obcore/math/linalg/eigen/Matrix.h"


namespace obvious
{

/**
 * @class SensorVelodyne3D
 * @brief class for velodyne 3D laser scanners (VLP16 PUCK, E32)
 * @author Jasmin Ziegler
 */
class SensorVelodyne3D : public Sensor
{
public:
  /**
   * Standard constructor
   * @param[in] raysIncl number of inclination rays of scanning device (vertical)
   * @param[in] inclMin lowest inclination angle in rad
   * @param[in] inclRes resolution of inclination rays in rad, i.e. angle between two vertical rays
   * @param[in] azimRes resolution of azimuth rays in rad, angle between two horizontal rays in 360° plane
   */
  SensorVelodyne3D(unsigned int raysIncl, double inclMin, double inclRes, double azimRes, double maxRange=INFINITY, double minRange=0.0, double lowReflectivityRange=INFINITY);

  /**
   * Destructor
   */
  virtual ~SensorVelodyne3D();

  /**
   * returns azimuth angle and inclination angle
   * @param[in] xCoord x coordinate of a 3D point
   * @param[in] yCoord y coordinate of a 3D point
   * @param[in] zCoord z coordinate of a 3D point
   * @param[out] inclAngle inclination angle in x-z-plane, 16 layers of rays from -15° to +15°, 2° resolution (VLP16 PUCK)
   * @param[out] azimAngle azimuth angle in x-y-plane, 0° to 360°
   */
  void returnAngles(double xCoord, double yCoord, double zCoord, double* inclAngle, double* azimAngle);

  /**
   * returns ray index
   * @param[in] azimAngle azimuth angle calculated in returnAngles()
   * @param[in] inclAngle inclination angle calculated in returnAngles()
   * @param[out] azimIndex index of azimuth ray number which is closest to 3D point to assign laser data to 3D point
   * @param[out] inclIndex index of inclination ray number which is closest to 3D point to assign laser data to 3D point
   * @todo make this function abstract so each velodyne sensor has to implement it since it may differ from sensor to sensor
   */
  void returnRayIndex(double azimAngle, double inclAngle, unsigned int* azimIndex, unsigned int* inclIndex);

  /**
   *
   */
  void setIndexMap(unsigned int raysAzim, unsigned int raysIncl);

  /**
   *
   */
  unsigned int lookupIndex(int indexSensormodel);

  /**
   * Project coordinate back to sensor index
   * @param[in] M matrix of 3D coordinates (homogeneous)
   * @param[out] indices vector of projection results (must be allocated outside)
   * @param[in] T temporary transformation matrix of coordinates
   */
  void backProject(obvious::Matrix* M, int* indices, obvious::Matrix* T=NULL);

private:

  double _azimRes;
  double _inclRes;
  int** _indexMap;
  /////////////////////////////////////////// DAS HIER SPÄTER RAUS UND IN AUFRUFENDER METHODE RICHTIG DIMENSIONIEREN!!!! getrows oder so
  int* _indices;
  //////////////////////////////


};

} /* namespace obvious */

#endif /* OBVISION_RECONSTRUCT_SPACE_SENSORVELODYNE3D_H_ */
