# Streamlit Community Cloud 部署说明

这个目录已经整理成适合上传到 GitHub 并部署到 Streamlit Community Cloud 的版本。

## 推荐仓库名
`cfrna-brain-tracing-app`

## 入口文件
`streamlit_app.py`

## 本地启动
```bash
pip install -r requirements.txt
streamlit run streamlit_app.py
```

## 部署到 Streamlit Community Cloud
1. 把整个目录推送到 GitHub 仓库。
2. 登录 `share.streamlit.io`。
3. 点击 **Create app**。
4. 选择你的 GitHub 仓库、分支和入口文件 `streamlit_app.py`。
5. 在 **App URL** 里可以填写你想要的子域名，例如：
   - `cfrna-brain-tracing`
   - `macaque-cfrna-tracing`
6. 点击部署。

## 推荐网址
可尝试填写：
`https://cfrna-brain-tracing.streamlit.app`

注意：这个网址只有在你自己的 Streamlit Community Cloud 账号中完成部署后才会真正生成。
如果该子域名已被占用，需要换一个。

## 当前包含内容
- Streamlit 前端
- SQLite 数据库文件 `cfrna_source_tracing.db`
- Bo2023 atlas 浏览器和重建矩阵接入口
- CLI 与 benchmark 代码

## 部署前建议
- 首次先在本地运行一次，确认页面能正常打开
- 如果后续要导入更大的重建矩阵，建议通过 GitHub 提交代码，不要直接在 Cloud 里手动修改
