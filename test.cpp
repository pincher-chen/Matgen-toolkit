#include <iostream>
#include <regex>
#include <string>
#include <fstream>
#include <regex>
#include <map>
#include "include/rapidjson/document.h"
#include "include/rapidjson/writer.h"
#include "include/rapidjson/stringbuffer.h"

using namespace rapidjson;
using namespace std;

struct Radius
{
    string species;
    string radius;
    string O_metal;

    Radius(string _species, string _radius, string _O_metal) : species(_species), radius(_radius), O_metal(_O_metal) {}

    bool operator<(const Radius &other) const
    {
        if (this->species < other.species)
        {
            return true;
        }
        return false;
    }

    Radius(const Radius &radius) {
        species = radius.species;
        O_metal = radius.O_metal;
        this->radius = radius.radius;
    }
};

class test {
public:
    void func() {
        Radius a = Radius("1", "1", "1");
        radius_dict.insert(pair<string, Radius>("1", a));
    }

    static map<string, Radius> radius_dict;
};

map<string, Radius> test::radius_dict;


int main()
{
    test t;
    t.func();
    // for(auto it = t.radius_dict.begin(); it != t.radius_dict.end(); it++) {
    //     cerr << it->first << endl;
    //     cerr << it->second.radius << endl;
    //     cerr << it->second.species << endl;
    //     cerr << it->second.O_metal << endl;
    // }
    auto it = t.radius_dict.find("1");
    cerr << it->first << endl;
    cerr << it->second.radius << endl;
    cerr << it->second.species << endl;
    cerr << it->second.O_metal << endl;

    //一，转json格式
    //1,获取Document对象
    Document doc;
    doc.SetObject(); //key-value 相当与map
    //doc.Setvalue();        //数组型 相当与vector
    Document::AllocatorType &allocator = doc.GetAllocator(); //获取分配器

    Value key;
    key.SetDouble(2.25);
    //2，给doc对象赋值
    doc.AddMember("name", key, allocator);
    Value& name = doc["name"];

    cout << name.GetDouble() << endl;
    //添加数组型数据
    // Value array1(kArrayType);
    // for (int i = 0; i < 3; i++)
    // {
    //     Value int_object(kObjectType);
    //     int_object.SetInt(i);
    //     array1.PushBack(int_object, allocator);
    // }

    // doc.AddMember("number", array1, allocator);

    // //3，将doc对象的值写入字符串
    // StringBuffer buffer;
    // //PrettyWriter<StringBuffer> writer(buffer);  //PrettyWriter是格式化的json，如果是Writer则是换行空格压缩后的json
    // Writer<StringBuffer> writer(buffer);
    // doc.Accept(writer);

    // cout << buffer.GetString() << endl;

    // //二，解析json格式
    // //1，将json格式字符串转换
    // string readdate;
    // readdate = buffer.GetString();
    // Document document;
    // document.Parse<0>(readdate.c_str());

    // //2,取出自己想要的值
    // Value &node1 = document["name"];
    // cout << "name:" << node1.GetString() << endl;

    // Value &node2 = document["number"];
    // cout << "number: " << endl;
    // if (node2.IsArray())
    // {
    //     for (int i = 0; i < node2.Size(); i++)
    //         cout << '\t' << node2[i].GetInt() << endl;
    // }

    return 0;
}