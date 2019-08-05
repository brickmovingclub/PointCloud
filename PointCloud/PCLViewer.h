#pragma once
class PCLViewer
{
public:
	PCLViewer();
	~PCLViewer();

	void ReadPcdFile(pcl::PointCloud<pcl::PointXYZRGBA>::Ptr cloud, const fs::path &fileName);		//	读取pcd文件中得点
	void AscToPcd(const fs::path &fileName);														//	将asc文件转成pcd文件格式
private:
	void ReadBuffer(const fs::path &fileDir, char **buffer);
private:
	unsigned int red;
	unsigned int green;
	unsigned int blue;
};

